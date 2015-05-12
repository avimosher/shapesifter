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
        std::pair<int,T> memory=force_memory[information->constraints[i]];
        information->value[i]=memory.second;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<const STORED_VOLUME_EXCLUSION_CONSTRAINT<T>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        force_memory[information->constraints[i]]=std::pair<int,T>(call_count,information->value[i]);
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment)
{
    call_count+=increment;
    LOG::cout<<"Count after increment: "<<call_count<<std::endl;
    auto information=std::static_pointer_cast<STORED_VOLUME_EXCLUSION_CONSTRAINT<T>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        if(force_memory.count(information->constraints[i])){// this could be compressed if I could be sure that the force would be initialized properly
            auto& memory=force_memory[information->constraints[i]];
            memory.first=call_count;
            memory.second+=increment*information->value[i];
            LOG::cout<<i<<" after increment: "<<memory.second<<std::endl;
        }
        else{
            force_memory[information->constraints[i]]=std::pair<int,T>(call_count,increment*information->value[i]);
        }
    }
    for(int i=0;i<constant_forces.size();i++){
        auto& memory=force_memory[constant_forces[i]];
        memory.first=call_count;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<T> rhs;
    RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
    constraints.clear();
    constant_forces.clear();
    if(stochastic){
        for(auto memory : force_memory){memory.second.second=(T)0;}}

    // let's say I rebuild the hierarchy every solve step
    // pass in list of indices and pre-built volumes    
    std::vector<AlignedBox<T,d>> bounding_list(rigid_data->Size());
    std::vector<int> index_list(rigid_data->Size());
    for(int s=0;s<rigid_data->Size();s++){
        auto structure=rigid_data->structures[s];
        bounding_list[s]=bounding_box(structure);
        index_list[s]=s;
    }
    KdBVH<T,3,int> hierarchy(index_list.begin(),index_list.end(),bounding_list.begin(),bounding_list.end());
    
    for(int s1=0;s1<rigid_data->structures.size();s1++){
        auto structure1=rigid_data->structures[s1];
        TEST_INTERSECTOR<TV> intersector(bounding_box(structure1));
        BVIntersect(hierarchy,intersector);
        //for(int s2=s1+1;s2<rigid_data->structures.size();s2++){
        for(int s2 : intersector.candidates){
            if(s1==s2){continue;}
            auto structure2=rigid_data->structures[s2];
            TV offset1,offset2;
            TV direction=structure1->Displacement(data,*structure2,offset1,offset2);
            T distance=direction.norm();
            T constraint_violation=distance-structure1->collision_radius-structure2->collision_radius;
            T slack_distance=-.005;
            T distance_condition=-.001;
            std::pair<int,T>& remembered=force_memory[CONSTRAINT(s1,s2)];
            LOG::cout<<"Constraint between "<<s1<<" and "<<s2<<" remembered force "<<remembered.second<<" count "<<remembered.first<<" "<<call_count<<" violation "<<constraint_violation<<" direction "<<direction.transpose()<<std::endl;
            if(constraint_violation<0){
                if(s1==0){
                    LOG::cout<<"Collision with "<<structure2->name<<" and "<<structure1->name<<std::endl;
                }
                TV x1=structure1->frame.position+offset1;
                TV x2=structure2->frame.position+offset2;
                TV object_offset1=structure1->frame.orientation.inverse()*offset1;
                TV object_offset2=structure2->frame.orientation.inverse()*offset2;
                LOG::cout<<"Offset1: "<<offset1<<" offset2: "<<offset2<<" x2-x1: "<<(x2-x1).transpose()<<std::endl;
                if(remembered.first!=call_count){
                    remembered.second=0;
                }
                CONSTRAINT_VECTOR dC_dA2=RIGID_STRUCTURE_INDEX_MAP<TV>::dC_dA(*structure2,object_offset2,x1,x2,direction);
                CONSTRAINT_VECTOR dC_dA1=RIGID_STRUCTURE_INDEX_MAP<TV>::dC_dA(*structure1,object_offset1,x1,x2,direction);
                T right_hand_side_force;
                if(constraint_violation<slack_distance){
                    LOG::cout<<"CONSTRAINT IS ON"<<std::endl;
                    right_hand_side_force=remembered.second;
                    terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s2,dC_dA2));
                    terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s1,-dC_dA1));
                    rhs.push_back(-constraint_violation+slack_distance);
#if 0
                    Matrix<T,t+d,t+d> force_balance_contribution2=remembered.second*RIGID_STRUCTURE_INDEX_MAP<TV>::dF_dA(*structure2,object_offset2,x1,x2,direction);
                    Matrix<T,t+d,t+d> force_balance_contribution1=remembered.second*RIGID_STRUCTURE_INDEX_MAP<TV>::dF_dA(*structure1,object_offset1,x1,x2,direction);
                
                    for(int j=0;j<t+d;j++){
                        for(int k=0;k<t+d;k++){
                            if(fabs(force_balance_contribution1(j,k))>1e-6){
                                force_terms.push_back(Triplet<T>(s1*(t+d)+j,s1*(t+d)+k,force_balance_contribution1(j,k)));
                            }
                            if(fabs(force_balance_contribution2(j,k))>1e-6){
                                force_terms.push_back(Triplet<T>(s2*(t+d)+j,s2*(t+d)+k,force_balance_contribution2(j,k)));
                            }
                        }
                    }
#endif
                    constraints.push_back(CONSTRAINT(s1,s2));
                }
                //else if(remembered.first==call_count && remembered.second>0){ // exponential falloff
                else if(remembered.first==call_count){ // exponential falloff
                    T exponent=-1/(1-sqr(constraint_violation/slack_distance-1));
                    right_hand_side_force=remembered.second*std::exp(1+exponent);
                    LOG::cout<<"SOFT CONSTRAINT IS ON"<<std::endl;
                    LOG::cout<<"EXP: "<<exponent<<" RHS contribution: "<<right_hand_side_force<<std::endl;
                    constant_forces.push_back(CONSTRAINT(s1,s2));

                    T constant_part=-2*right_hand_side_force*sqr(exponent)*(constraint_violation/slack_distance-1)/slack_distance;
                    auto dc_dx1=RIGID_STRUCTURE_INDEX_MAP<TV>::dXN_dA(*structure1,object_offset1,x1,x2);
                    auto dc_dx2=RIGID_STRUCTURE_INDEX_MAP<TV>::dXN_dA(*structure2,object_offset2,x1,x2);
                    Matrix<T,t+d,t+d> dA1dx1=dC_dA1.transpose()*constant_part*dc_dx1;
                    Matrix<T,t+d,t+d> dA1dx2=-dC_dA1.transpose()*constant_part*dc_dx2;
                    Matrix<T,t+d,t+d> dA2dx1=-dC_dA2.transpose()*constant_part*dc_dx1;
                    Matrix<T,t+d,t+d> dA2dx2=dC_dA2.transpose()*constant_part*dc_dx2;
                    Flatten_Matrix_Term(s1,s1,dA1dx1,force_terms);
                    Flatten_Matrix_Term(s1,s2,dA1dx2,force_terms);
                    Flatten_Matrix_Term(s2,s1,dA2dx1,force_terms);
                    Flatten_Matrix_Term(s2,s2,dA2dx2,force_terms);
                }
                right_hand_side.template block<t+d,1>(s1*(t+d),0)+=dC_dA1.transpose()*right_hand_side_force;
                right_hand_side.template block<t+d,1>(s2*(t+d),0)-=dC_dA2.transpose()*right_hand_side_force;
            }
        }
    }
    constraint_terms.resize(constraints.size(),rigid_data->Velocity_DOF());
    constraint_rhs.resize(rhs.size(),1);
    for(int i=0;i<rhs.size();i++){constraint_rhs(i,0)=rhs[i];}
    Flatten_Matrix(terms,constraint_terms);
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
    LOG::cout<<"Volume exclusion constraints: "<<constraints.size()<<std::endl;
    for(int i=0;i<constraints.size();i++){
        auto lineGeometry=new osg::Geometry();
        auto vertices=new osg::Vec3Array(2);
        const CONSTRAINT& constraint=constraints[i];
        int body_index1=constraint.first;
        int body_index2=constraint.second;
        LOG::cout<<"Between "<<body_index1<<" and "<<body_index2<<std::endl;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        auto firstAttachment=rigid_structure1->frame.position;
        auto secondAttachment=rigid_structure2->frame.position;
        (*vertices)[0].set(firstAttachment(0),firstAttachment(1),firstAttachment(2));
        (*vertices)[1].set(secondAttachment(0),secondAttachment(1),secondAttachment(2));
        lineGeometry->setVertexArray(vertices);
        auto colors=new osg::Vec4Array;
        colors->push_back(osg::Vec4(0.0f,0.0f,1.0f,1.0f));
        lineGeometry->setColorArray(colors);
        lineGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
        auto normals=new osg::Vec3Array;
        normals->push_back(osg::Vec3f(0.0f,-1.0f,0.0f));
        lineGeometry->setNormalArray(normals);
        lineGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);
        
        lineGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,2));
        auto lineGeode=new osg::Geode();
        lineGeode->addDrawable(lineGeometry);
        volume_exclusion_group->addChild(lineGeode);
    }
    group->addChild(volume_exclusion_group);
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(VOLUME_EXCLUSION_CONSTRAINT)
GENERIC_TYPE_DEFINITION(VOLUME_EXCLUSION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(VOLUME_EXCLUSION_CONSTRAINT,void)
{
    auto volume_exclusion_constraint=simulation.force.template Find_Or_Create<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
    return 0;
}
