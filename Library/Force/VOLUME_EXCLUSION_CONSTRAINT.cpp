#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/VOLUME_EXCLUSION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/HASHING.h>
#include <Utilities/LOG.h>
#include <Utilities/OSG_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
#include <unsupported/Eigen/BVH>
#include <osg/Geometry>
#include <osg/Geode>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> std::shared_ptr<FORCE_REFERENCE<typename TV::Scalar>> VOLUME_EXCLUSION_CONSTRAINT<TV>::
Create_Stored_Force() const
{
    return std::static_pointer_cast<FORCE_REFERENCE<T>>(std::make_shared<STORED_VOLUME_EXCLUSION_CONSTRAINT<T>>());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<STORED_VOLUME_EXCLUSION_CONSTRAINT<T>>(force_information);
    information->constraints=constraints;
    information->value.resize(constraints.size());
    for(int i=0;i<information->constraints.size();i++){
        information->value[i]=force_memory[information->constraints[i]].second;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<const STORED_VOLUME_EXCLUSION_CONSTRAINT<T>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        force_memory[information->constraints[i]]=std::pair<int,T>(call_count,information->value[i]);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment)
{
    call_count+=increment;
    auto information=std::static_pointer_cast<STORED_VOLUME_EXCLUSION_CONSTRAINT<T>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        if(force_memory.count(information->constraints[i])){// this could be compressed if I could be sure that the force would be initialized properly
            auto& memory=force_memory[information->constraints[i]];
            memory.first=call_count;
            memory.second+=increment*information->value[i];}
        else{
            force_memory[information->constraints[i]]=std::pair<int,T>(call_count,increment*information->value[i]);}
    }
    for(int i=0;i<constant_forces.size();i++){
        std::get<0>(constant_force_memory[constant_forces[i]])=call_count;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    typedef Matrix<T,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE,1> FORCE_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<Triplet<FORCE_VECTOR>> forces;
    std::vector<T> rhs;
    constraints.clear();
    constant_forces.clear();
    if(stochastic){
        for(auto& memory : force_memory){memory.second.second=(T)0;memory.second.first=-1;}
        for(auto& constant_memory : constant_force_memory){std::get<2>(constant_memory.second)=(T)0;std::get<0>(constant_memory.second)=-1;}}

    std::vector<AlignedBox<T,d>> bounding_list(rigid_data->Size());
    std::vector<int> index_list(rigid_data->Size());
    for(int s=0;s<rigid_data->Size();s++){
        auto structure=rigid_data->structures[s];
        bounding_list[s]=structure->Bounding_Box();
        index_list[s]=s;}
    KdBVH<T,3,int> hierarchy(index_list.begin(),index_list.end(),bounding_list.begin(),bounding_list.end());
    
    for(int s1=0;s1<rigid_data->structures.size();s1++){
        auto structure1=rigid_data->structures[s1];
        BOX_PROXIMITY_SEARCH<TV> intersector(data,structure1->Bounding_Box());
        BVIntersect(hierarchy,intersector);
        for(int s2 : intersector.candidates){
            if(s1>=s2){continue;}
            auto structure2=rigid_data->structures[s2];
            std::vector<TV> offsets(2);
            TV relative_position=structure1->Displacement(data,*structure2,offsets[0],offsets[1]);
            TV direction=relative_position.normalized();
            T constraint_violation=relative_position.norm()-structure1->collision_radius-structure2->collision_radius;
            //LOG::cout<<"Relative position: "<<relative_position.transpose()<<" normalized: "<<direction.transpose()<<" offset 1: "<<offsets[0].transpose()<<" offset 2: "<<offsets[1].transpose()<<std::endl;
            T slack_distance=-.005;
            T push_out_distance=1e-8;
            if(constraint_violation<0){
                CONSTRAINT constraint(s1,s2);
                std::pair<int,T>& memory=force_memory[constraint];
                auto& constant_memory=constant_force_memory[constraint];
                std::vector<T_SPIN> spins={structure1->twist.angular,structure2->twist.angular};
                CONSTRAINT_VECTOR dC_dA1=RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(spins[0],offsets[0],relative_position);
                CONSTRAINT_VECTOR dC_dA2=RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(spins[1],offsets[1],relative_position);
                T right_hand_side_force=0;
                FORCE_VECTOR force_direction1=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(offsets[0]).transpose()*direction;
                FORCE_VECTOR force_direction2=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(offsets[1]).transpose()*direction;
                if(constraint_violation<slack_distance){
                    if(std::get<0>(memory)!=call_count){
                        std::get<1>(memory)=0;
                        if(std::get<0>(constant_memory)==call_count){
                            std::get<1>(memory)=std::get<1>(constant_memory)*sqr(constraint_violation);}}
                    LOG::cout<<"Constraint between "<<s1<<" and "<<s2<<" with force: "<<std::get<1>(memory)<<std::endl;
                    right_hand_side_force=std::get<1>(memory);
                    terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s2,dC_dA2));
                    terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s1,-dC_dA1));
                    forces.push_back(Triplet<FORCE_VECTOR>(s2,constraints.size(),force_direction2));
                    forces.push_back(Triplet<FORCE_VECTOR>(s1,constraints.size(),-force_direction1));
                    rhs.push_back(-constraint_violation+slack_distance+push_out_distance);
                    constraints.push_back(constraint);}
                else if(memory.first==call_count || std::get<0>(constant_memory)==call_count){ // exponential falloff
                    if(std::get<0>(constant_memory)!=call_count){
                        std::get<1>(constant_memory)=std::get<1>(memory)/sqr(constraint_violation);}

                    right_hand_side_force=std::get<1>(constant_memory)*sqr(constraint_violation);
                    constant_forces.push_back(constraint);
                    //LOG::cout<<"Force between "<<s1<<" and "<<s2<<" with force: "<<std::get<2>(constant_memory)<<" but also "<<right_hand_side_force<<" characteristic length "<<characteristic_length<<std::endl;

                    T constant_part=-2*right_hand_side_force*sqr(exponent)*(constraint_violation/characteristic_length-1)/characteristic_length;
                    //T constant_part=right_hand_side_force*sqr(exponent)*(-2*(constraint_violation/characteristic_length-1))/characteristic_length
                    std::vector<int> indices={s1,s2};
                    for(int s1=0;s1<2;s1++){
                        for(int s2=0;s2<2;s2++){
                            int term_sign=s1==s2?1:-1;
                            Flatten_Matrix_Term<T,t+d,t+d,d,d>(indices[s1],indices[s2],0,0,RIGID_STRUCTURE_INDEX_MAP<TV>::dForce_dVelocity(relative_position)*right_hand_side_force*term_sign,force_terms);
                            Flatten_Matrix_Term<T,t+d,t+d,d,t>(indices[s1],indices[s2],0,1,RIGID_STRUCTURE_INDEX_MAP<TV>::dForce_dSpin(relative_position,spins[s2],offsets[s2])*right_hand_side_force*term_sign,force_terms);
                            Flatten_Matrix_Term<T,t+d,t+d,t,d>(indices[s1],indices[s2],1,0,RIGID_STRUCTURE_INDEX_MAP<TV>::dTorque_dVelocity(relative_position,offsets[s1])*right_hand_side_force*term_sign,force_terms);
                            Flatten_Matrix_Term<T,t+d,t+d,t,t>(indices[s1],indices[s2],1,1,RIGID_STRUCTURE_INDEX_MAP<TV>::dTorque_dSpin(relative_position,s1,s2,spins[s2],offsets)*right_hand_side_force*term_sign,force_terms);
                        }}

                    // TODO: these derivatives are probably wrong.
#if 0
                    Matrix<T,t+d,t+d> dA1dx1=dC_dA1.transpose()*constant_part*dC_dA1;
                    Matrix<T,t+d,t+d> dA1dx2=-dC_dA1.transpose()*constant_part*dC_dA2;
                    Matrix<T,t+d,t+d> dA2dx1=-dC_dA2.transpose()*constant_part*dC_dA1;
                    Matrix<T,t+d,t+d> dA2dx2=dC_dA2.transpose()*constant_part*dC_dA2;
                    Flatten_Matrix_Term(s1,s1,dA1dx1,force_terms);
                    Flatten_Matrix_Term(s1,s2,dA1dx2,force_terms);
                    Flatten_Matrix_Term(s2,s1,dA2dx1,force_terms);
                    Flatten_Matrix_Term(s2,s2,dA2dx2,force_terms);
#endif
                }
                LOG::cout<<"Force applied to body "<<s1<<": "<<(force_direction1*right_hand_side_force).transpose()<<std::endl;
                LOG::cout<<"Force applied to body "<<s2<<": "<<(-force_direction2*right_hand_side_force).transpose()<<std::endl;
                right_hand_side.template block<t+d,1>(s1*(t+d),0)+=force_direction1*right_hand_side_force;
                right_hand_side.template block<t+d,1>(s2*(t+d),0)-=force_direction2*right_hand_side_force;}}}
    constraint_rhs.resize(rhs.size(),1);
    for(int i=0;i<rhs.size();i++){constraint_rhs(i,0)=rhs[i];}
    constraint_terms.resize(constraints.size(),rigid_data->Velocity_DOF());
    Flatten_Matrix(terms,constraint_terms);
    constraint_forces.resize(rigid_data->Velocity_DOF(),constraints.size());
    Flatten_Matrix(forces,constraint_forces);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
    osg::Group* group=node->asGroup();
    group->removeChild(getNamedChild(group,Static_Name()));
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    osg::Group* volume_exclusion_group=new osg::Group();
    volume_exclusion_group->setName(Static_Name());
    for(int i=0;i<constraints.size();i++){
        auto lineGeometry=new osg::Geometry();
        auto vertices=new osg::Vec3Array(2);
        const CONSTRAINT& constraint=constraints[i];
        int body_index1=constraint.first;
        int body_index2=constraint.second;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        auto firstAttachment=rigid_structure1->frame.position;
        auto secondAttachment=rigid_structure2->frame.position;
        (*vertices)[0].set(firstAttachment(0),firstAttachment(1),firstAttachment(2));
        (*vertices)[1].set(secondAttachment(0),secondAttachment(1),secondAttachment(2));
        lineGeometry->setVertexArray(vertices);
        auto colors=new osg::Vec4Array;
        colors->push_back(osg::Vec4(1.0f,0.0f,0.0f,1.0f));
        lineGeometry->setColorArray(colors);
        lineGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
        auto normals=new osg::Vec3Array;
        normals->push_back(osg::Vec3f(0.0f,-1.0f,0.0f));
        lineGeometry->setNormalArray(normals);
        lineGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);
        
        lineGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,2));

        osg::StateSet* stateset=new osg::StateSet;
        stateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
        lineGeometry->setStateSet(stateset);

        auto lineGeode=new osg::Geode();
        lineGeode->addDrawable(lineGeometry);
        volume_exclusion_group->addChild(lineGeode);}
    group->addChild(volume_exclusion_group);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(VOLUME_EXCLUSION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(VOLUME_EXCLUSION_CONSTRAINT,void)
{
    auto volume_exclusion_constraint=simulation.force.template Find_Or_Create<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
    return 0;
}
