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
#ifdef VIEWER
#include <osg/Geometry>
#include <osg/Geode>
#endif
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
    information->constraints.resize(constraint_force_indices.size());
    information->value.resize(constraint_force_indices.size());
    for(int i=0;i<information->constraints.size();i++){
        information->constraints[i]=constraint_indices[constraint_force_indices[i]];
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
        std::get<MEMORY_COUNT>(constant_force_memory[constant_force_indices[i]])=call_count;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Identify_Interactions_And_Compute_Errors(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto& structures=rigid_data->structures;
    Matrix<T,Dynamic,1>& right_hand_side=system.RHS(data,force,*rigid_data);
    Matrix<T,Dynamic,1>& constraint_right_hand_side=system.RHS(data,force,*this);
    int new_constraints=0,old_constraints=0;
    constraints.clear();constant_forces.clear();constraint_indices.clear();constant_force_indices.clear();rhs.clear();
    constraint_force_indices.clear();
    if(stochastic){
        constraint_count.push(0);
        // TODO: keep a list of which memories were touched and clean only those
        for(auto& memory : force_memory){std::get<MEMORY_FORCE>(memory.second)=(T)0;std::get<MEMORY_COUNT>(memory.second)=-1;}
        for(auto& constant_memory : constant_force_memory){std::get<MEMORY_FORCE>(constant_memory.second)=(T)0;std::get<MEMORY_COUNT>(constant_memory.second)=-1;}}

    int substructure_count=0;
    for(int s=0;s<rigid_data->Size();s++){
        substructure_count+=structures[s]->Substructures().size();}
    std::vector<AlignedBox<T,d>> bounding_list(substructure_count);
    std::vector<std::pair<int,int>> index_list(substructure_count);
    int current_substructure=0;
    for(int s=0;s<rigid_data->Size();s++){
        auto structure=structures[s];
        for(int i=0;i<structure->Substructures().size();i++,current_substructure++){
            bounding_list[current_substructure]=structure->Substructure(i).Bounding_Box(structure->frame);
            index_list[current_substructure]=std::pair<int,int>(s,i);}
    }
    KdBVH<T,3,std::pair<int,int>> hierarchy(index_list.begin(),index_list.end(),bounding_list.begin(),bounding_list.end());
    
    std::array<TV,2> offsets;
    std::array<FORCE_VECTOR,2> force_directions;
    LOG::cout<<"Before loop"<<std::endl;
    for(int s1=0;s1<structures.size();s1++){
        auto structure1=structures[s1];
        for(int ss1=0;ss1<structure1->Substructures().size();ss1++){
            auto& substructure1=structure1->Substructure(ss1);
            BOX_PROXIMITY_SEARCH<TV> intersector(data,substructure1.Bounding_Box(structure1->frame));
            BVIntersect(hierarchy,intersector);
            for(const std::pair<int,int>& s2 : intersector.candidates){
                if(s1>=s2.first){continue;}
                auto structure2=structures[s2.first];
                int ss2=s2.second;
                auto& substructure2=structure2->Substructure(ss2);
                TV relative_position=SUBSTRUCTURE<TV>::Displacement(data,structure1->frame,structure2->frame,substructure1,substructure2,offsets);
                T threshold=substructure1.radius+substructure2.radius;
                T constraint_violation=relative_position.norm()-threshold;
                T slack_distance=-.005,push_out_distance=1e-8;
                LOG::cout<<"s1: "<<s1<<" ss1: "<<ss1<<" s2: "<<s2.first<<" ss2: "<<ss2<<" viol: "<<constraint_violation<<std::endl;
                if(constraint_violation<0){
                    INDICES indices={s1,s2.first,ss1,ss2};
                    auto& memory=force_memory[indices];
                    auto& constant_memory=constant_force_memory[indices];
                    std::array<T_SPIN,2> spins={structure1->twist.angular,structure2->twist.angular};
                    T right_hand_side_force=0;
                    for(int s=0;s<2;s++){force_directions[s]=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(offsets[s]).transpose()*relative_position.normalized();}
                    if(constraint_violation<slack_distance){
                        if(std::get<MEMORY_COUNT>(memory)!=call_count){ // if the force wasn't on last step, start it at zero
                            new_constraints++;
                            std::get<MEMORY_FORCE>(memory)=0;
                            if(std::get<MEMORY_COUNT>(constant_memory)==call_count){ // if we're moving from a penalty force, use its magnitude
                                std::get<MEMORY_FORCE>(memory)=std::get<MEMORY_FORCE>(constant_memory)*sqr(constraint_violation);}}
                        else{old_constraints++;}
                        right_hand_side_force=std::get<MEMORY_FORCE>(memory);
                        /*if(right_hand_side_force<0 && constraint_violation<remembered_threshold){
                            // TODO: switch the constraint force back on
                            }*/
                        if(right_hand_side_force>=0){
                            constraint_force_indices.push_back(constraint_indices.size());
                        }
                        rhs.push_back(constraint_violation-slack_distance-push_out_distance);
                        constraints.push_back(std::make_tuple(force_directions,spins,offsets,relative_position));
                        constraint_indices.push_back(indices);}
                    if((constraint_violation>=slack_distance || right_hand_side_force<0) && (std::get<MEMORY_COUNT>(constant_memory)==call_count || std::get<MEMORY_COUNT>(memory)==call_count)){
                        if(std::get<MEMORY_COUNT>(memory)==call_count){ // this will hold whether the remembered force is positive or negative
                                std::get<MEMORY_FORCE>(constant_memory)=std::get<MEMORY_FORCE>(memory)/sqr(constraint_violation);
                        }
                        //if(std::get<MEMORY_COUNT>(constant_memory)!=call_count){std::get<MEMORY_FORCE>(constant_memory)=std::get<MEMORY_FORCE>(memory)/sqr(constraint_violation);}
                        right_hand_side_force=std::get<MEMORY_FORCE>(constant_memory)*sqr(constraint_violation);
                        CONSTANT_FORCE constant_force(spins,offsets,relative_position,threshold);
                        constant_forces.push_back(constant_force);
                        constant_force_indices.push_back(indices);}
                    LOG::cout<<"RHS force: "<<right_hand_side_force<<std::endl;
                    right_hand_side.template block<t+d,1>(s1*(t+d),0)-=force_directions[0]*right_hand_side_force;
                    right_hand_side.template block<t+d,1>(s2.first*(t+d),0)+=force_directions[1]*right_hand_side_force;}}}}
    equations_changed=new_constraints>0 || old_constraints!=constraint_count.peek();
    constraint_right_hand_side=Map<Matrix<T,Dynamic,1>>(rhs.data(),rhs.size());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Compute_Derivatives(DATA<TV>& data,FORCE<TV>& force,MATRIX_BUNDLE<TV>& system)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    std::vector<Triplet<T>>& force_terms=system.Matrix_Block_Terms(data,force,*rigid_data);
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<Triplet<FORCE_VECTOR>> forces;

    for(int i=0;i<constraints.size();i++){
        auto& constraint=constraints[i];
        auto& indices=constraint_indices[i];
        auto& memory=force_memory[indices];
        for(int s=0,sgn=-1;s<2;s++,sgn+=2){
            terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,indices[s],sgn*RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(std::get<CONSTRAINT_SPINS>(constraint)[s],std::get<CONSTRAINT_OFFSETS>(constraint)[s],std::get<CONSTRAINT_RELATIVE_POSITION>(constraint))));
            /*forces.push_back(Triplet<FORCE_VECTOR>(indices[s],i,sgn*std::get<CONSTRAINT_FORCE_DIRECTIONS>(constraint)[s]));*/}
        RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Constraint_Force_Derivatives({indices[0],indices[1]},std::get<MEMORY_FORCE>(memory),std::get<CONSTRAINT_RELATIVE_POSITION>(constraint),std::get<CONSTRAINT_OFFSETS>(constraint),std::get<CONSTRAINT_SPINS>(constraint),force_terms);}

    for(int i=0;i<constraint_force_indices.size();i++){
        auto& constraint=constraints[constraint_force_indices[i]];
        auto& indices=constraint_indices[constraint_force_indices[i]];
        for(int s=0,sgn=-1;s<2;s++,sgn+=2){
            forces.push_back(Triplet<FORCE_VECTOR>(indices[s],i,sgn*std::get<CONSTRAINT_FORCE_DIRECTIONS>(constraint)[s]));}
    }

    for(int i=0;i<constant_forces.size();i++){
        auto& constant_force=constant_forces[i];
        auto& indices=constant_force_indices[i];
        auto& constant_memory=constant_force_memory[indices];
        RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Penalty_Force_Derivatives({indices[0],indices[1]},std::get<CONSTANT_FORCE_THRESHOLD>(constant_force),std::get<MEMORY_FORCE>(constant_memory),std::get<CONSTANT_FORCE_RELATIVE_POSITION>(constant_force),std::get<CONSTANT_FORCE_OFFSETS>(constant_force),std::get<CONSTANT_FORCE_SPINS>(constant_force),force_terms);}

    system.Flatten_Jacobian_Block(data,force,*this,*rigid_data,terms);
    system.Flatten_Jacobian_Block(data,force,*rigid_data,*this,forces);
}
///////////////////////////////////////////////////////////////////////
#ifdef VIEWER
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
        const auto& constraint=constraint_indices[i];
        auto rigid_structure1=rigid_data->structures[constraint[0]];
        auto rigid_structure2=rigid_data->structures[constraint[1]];
        auto firstAttachment=rigid_structure1->frame.position;
        auto secondAttachment=rigid_structure2->frame.position;
        osg::Vec4 color(1.0f,0.0f,0.0f,1.0f);
        auto lineGeode=createLine(color);
        updateLine<TV>(lineGeode,{firstAttachment,secondAttachment});
        volume_exclusion_group->addChild(lineGeode);

        if(errors.rows()==constraints.size()){
            osg::Vec4 errorColor(1.0f,1.0f,1.0f,1.0f);
            auto errorLineGeode=createLine(errorColor);
            TV mean=(firstAttachment+secondAttachment)/2;
            TV offset=mean+TV::Unit(2)*errors(i);
            updateLine<TV>(errorLineGeode,{mean,offset});
            volume_exclusion_group->addChild(errorLineGeode);}
    }
    group->addChild(volume_exclusion_group);
}
#endif
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(VOLUME_EXCLUSION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(VOLUME_EXCLUSION_CONSTRAINT,void)
{
    auto volume_exclusion_constraint=simulation.force.template Find_Or_Create<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
    return 0;
}
