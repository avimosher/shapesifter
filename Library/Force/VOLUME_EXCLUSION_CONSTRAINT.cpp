#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/MATRIX_BUNDLE.h>
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
    if(increment==1){
        constraint_count.push(information->constraints.size());}
    else{
        constraint_count.pop();}
    for(int i=0;i<information->constraints.size();i++){
        if(force_memory.count(information->constraints[i])){// this could be compressed if I could be sure that the force would be initialized properly
            auto& memory=force_memory[information->constraints[i]];
            memory.first=call_count;
            memory.second+=increment*information->value[i];}
        else{
            force_memory[information->constraints[i]]=std::pair<int,T>(call_count,increment*information->value[i]);}}
    for(int i=0;i<constant_forces.size();i++){
        std::get<0>(constant_force_memory[constant_forces[i]])=call_count;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    Matrix<T,Dynamic,1>& right_hand_side=system.RHS(data,force,*rigid_data);
    Matrix<T,Dynamic,1>& constraint_right_hand_side=system.RHS(data,force,*this);
    std::vector<Triplet<T>>& force_terms=system.Matrix_Block_Terms(data,force,*rigid_data);
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<Triplet<FORCE_VECTOR>> forces;
    int new_constraints=0;
    int old_constraints=0;
    constraints.clear();
    constant_forces.clear();
    rhs.clear();
    if(stochastic){
        constraint_count.push(0);
        for(auto& memory : force_memory){std::get<1>(memory.second)=(T)0;std::get<0>(memory.second)=-1;}
        for(auto& constant_memory : constant_force_memory){std::get<1>(constant_memory.second)=(T)0;std::get<0>(constant_memory.second)=-1;}}

    std::vector<AlignedBox<T,d>> bounding_list(rigid_data->Size());
    std::vector<int> index_list(rigid_data->Size());
    for(int s=0;s<rigid_data->Size();s++){
        auto structure=rigid_data->structures[s];
        bounding_list[s]=structure->Bounding_Box();
        index_list[s]=s;}
    KdBVH<T,3,int> hierarchy(index_list.begin(),index_list.end(),bounding_list.begin(),bounding_list.end());
    
    std::vector<TV> offsets(2);
    std::vector<FORCE_VECTOR> force_directions(2);
    for(int s1=0;s1<rigid_data->structures.size();s1++){
        auto structure1=rigid_data->structures[s1];
        BOX_PROXIMITY_SEARCH<TV> intersector(data,structure1->Bounding_Box());
        BVIntersect(hierarchy,intersector);
        for(int s2 : intersector.candidates){
            if(s1>=s2){continue;}
            auto structure2=rigid_data->structures[s2];
            TV relative_position=structure1->Displacement(data,*structure2,offsets[0],offsets[1]);
            TV direction=relative_position.normalized();
            T threshold=structure1->collision_radius+structure2->collision_radius;
            T constraint_violation=relative_position.norm()-threshold;
            T slack_distance=-.005;
            T push_out_distance=1e-8;
            if(constraint_violation<0){
                CONSTRAINT constraint(s1,s2);
                auto& memory=force_memory[constraint];
                auto& constant_memory=constant_force_memory[constraint];
                std::vector<T_SPIN> spins={structure1->twist.angular,structure2->twist.angular};
                T right_hand_side_force=0;
                for(int s=0;s<2;s++){force_directions[s]=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(offsets[s]).transpose()*direction;}
                std::vector<int> indices={s1,s2};
                if(constraint_violation<slack_distance){
                    if(std::get<0>(memory)!=call_count){ // if the force wasn't on last step, start it at zero
                        new_constraints++;
                        std::get<1>(memory)=0;
                        if(std::get<0>(constant_memory)==call_count){ // if we're moving from a penalty force, use its magnitude
                            std::get<1>(memory)=std::get<1>(constant_memory)*sqr(constraint_violation);}}
                    else{
                        old_constraints++;}
                    right_hand_side_force=std::get<1>(memory);
                    for(int s=0,sgn=-1;s<2;s++,sgn+=2){
                        terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),indices[s],sgn*RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(spins[s],offsets[s],relative_position)));
                        forces.push_back(Triplet<FORCE_VECTOR>(indices[s],constraints.size(),sgn*force_directions[s]));}
                    RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Constraint_Force_Derivatives(indices,right_hand_side_force,relative_position,offsets,spins,force_terms);
                    rhs.push_back(-constraint_violation+slack_distance+push_out_distance);
                    constraints.push_back(constraint);}
                else if(std::get<0>(memory)==call_count || std::get<0>(constant_memory)==call_count){
                    if(std::get<0>(constant_memory)!=call_count){
                        std::get<1>(constant_memory)=std::get<1>(memory)/sqr(constraint_violation);}
                    right_hand_side_force=std::get<1>(constant_memory)*sqr(constraint_violation);
                    constant_forces.push_back(constraint);
                    RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Penalty_Force_Derivatives(indices,threshold,std::get<1>(constant_memory),relative_position,offsets,spins,force_terms);
                }
                right_hand_side.template block<t+d,1>(s1*(t+d),0)+=force_directions[0]*right_hand_side_force;
                right_hand_side.template block<t+d,1>(s2*(t+d),0)-=force_directions[1]*right_hand_side_force;}}}
    LOG::cout<<"Existing constraints: "<<constraint_count.peek()<<" old: "<<old_constraints<<" new: "<<new_constraints<<std::endl;
    equations_changed=new_constraints>0 || old_constraints!=constraint_count.peek();
    constraint_right_hand_side.resize(rhs.size(),1);
    for(int i=0;i<rhs.size();i++){constraint_right_hand_side(i,0)=rhs[i];}

    system.Flatten_Jacobian_Block(data,force,*this,*rigid_data,terms);
    system.Flatten_Jacobian_Block(data,force,*rigid_data,*this,forces);
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
        std::vector<TV> points={firstAttachment,secondAttachment};
        auto lineGeode=createLine(osg::Vec4(1.0f,0.0f,0.0f,1.0f));
        updateLine(lineGeode,points);
        volume_exclusion_group->addChild(lineGeode);

        if(errors.rows()==constraints.size()){
            auto errorLineGeode=createLine(osg::Vec4(1.0f,1.0f,1.0f,1.0f));
            auto mean=(firstAttachment+secondAttachment)/2;
            auto offset=mean+TV::Unit(2)*errors(i);
            std::vector<TV> error_points={mean,offset};
            updateLine(errorLineGeode,error_points);
            volume_exclusion_group->addChild(errorLineGeode);}
    }
    group->addChild(volume_exclusion_group);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(VOLUME_EXCLUSION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(VOLUME_EXCLUSION_CONSTRAINT,void)
{
    auto volume_exclusion_constraint=simulation.force.template Find_Or_Create<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
    return 0;
}
