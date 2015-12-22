#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Equation/MATRIX_BUNDLE.h>
#include <Force/FORCE.h>
#include <Force/WALL_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> std::shared_ptr<FORCE_REFERENCE<typename TV::Scalar>> WALL_CONSTRAINT<TV>::
Create_Stored_Force() const
{
    return std::static_pointer_cast<FORCE_REFERENCE<T>>(std::make_shared<STORED_WALL_CONSTRAINT>());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void WALL_CONSTRAINT<TV>::
Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<STORED_WALL_CONSTRAINT>(force_information);
    information->constraints=constraints;
    information->value.resize(constraints.size());
    for(int i=0;i<information->constraints.size();i++){
        information->value[i]=force_memory[information->constraints[i]].second;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void WALL_CONSTRAINT<TV>::
Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<const STORED_WALL_CONSTRAINT>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        force_memory[information->constraints[i]]=std::pair<int,T>(call_count,information->value[i]);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void WALL_CONSTRAINT<TV>::
Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment)
{
    call_count+=increment;
    auto information=std::static_pointer_cast<STORED_WALL_CONSTRAINT>(force_information);
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
template<class TV> void WALL_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    if(stochastic){
        for(auto& memory : force_memory){std::get<1>(memory.second)=(T)0;std::get<0>(memory.second)=-1;}
        for(auto& constant_memory : constant_force_memory){std::get<1>(constant_memory.second)=(T)0;std::get<0>(constant_memory.second)=-1;}}

    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    Matrix<T,Dynamic,1>& right_hand_side=system.RHS(data,force,*rigid_data);
    Matrix<T,Dynamic,1>& constraint_right_hand_side=system.RHS(data,force,*this);
    std::vector<Triplet<T>>& force_terms=system.Matrix_Block_Terms(data,force,*rigid_data);
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<Triplet<FORCE_VECTOR>> forces;

    constraints.clear();
    constant_forces.clear();
    std::vector<T> rhs;
    T slack_distance=-.005;
    T push_out_distance=1e-8;
    for(int s=0;s<rigid_data->structures.size();s++){
        auto structure=rigid_data->structures[s];
        T_SPIN spin=structure->twist.angular;
        std::vector<TV> extrema={structure->frame.orientation*(-structure->collision_extent),structure->frame.orientation*structure->collision_extent};
        std::vector<TV> axis_extrema(2);
        T threshold=structure->collision_radius;
        // check against each wall
        for(int axis=0;axis<d;axis++){
            TV axis_minimum,axis_maximum;
            if(extrema[0][axis]>extrema[1][axis]){
                axis_extrema[0]=extrema[0];
                axis_extrema[1]=extrema[1];}
            else if (extrema[0][axis]<extrema[1][axis]){
                axis_extrema[0]=extrema[1];
                axis_extrema[1]=extrema[0];}
            else{
                axis_extrema[0]=axis_extrema[1]=(extrema[0]+extrema[1])*(T).5;}

            for(int w=0,sgn=1;w<2;w++,sgn-=2){
                if(walls[w][axis]){
                    T relative_position=sgn*(axis_extrema[w][axis]+structure->frame.position[axis]-data.domain[w][axis]);
                    T constraint_violation=relative_position-threshold;
                    TV offset=axis_extrema[w];
                    TV direction=sgn*TV::Unit(axis);
                    if(constraint_violation<0){ // do something about it
                        CONSTRAINT constraint(s,axis,w);
                        auto& memory=force_memory[constraint];
                        auto& constant_memory=constant_force_memory[constraint];
                        T right_hand_side_force=0;
                        FORCE_VECTOR force_direction=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(offset).transpose()*direction;
                        if(constraint_violation<slack_distance){
                            if(std::get<0>(memory)!=call_count){
                                std::get<1>(memory)=0;
                                if(std::get<0>(constant_memory)==call_count){
                                    std::get<1>(memory)=std::get<1>(constant_memory)*sqr(constraint_violation);}}
                            //LOG::cout<<"Constraint force with base: "<<std::get<1>(memory)<<std::endl;
                            right_hand_side_force=std::get<1>(memory);
                            terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s,RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(spin,offset,direction*relative_position,slack_distance+push_out_distance)));
                            forces.push_back(Triplet<FORCE_VECTOR>(s,constraints.size(),force_direction));
                            RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Constraint_Force_Derivative(s,right_hand_side_force,direction*relative_position,offset,spin,force_terms);
                            rhs.push_back(-constraint_violation+slack_distance+push_out_distance);
                            constraints.push_back(constraint);}
                        else if(std::get<0>(memory)==call_count || std::get<0>(constant_memory)==call_count){
                            if(std::get<0>(constant_memory)!=call_count){
                                std::get<1>(constant_memory)=std::get<1>(memory)/sqr(constraint_violation);}
                            //LOG::cout<<"Penalty force with base: "<<std::get<1>(constant_memory)<<std::endl;
                            right_hand_side_force=std::get<1>(constant_memory)*sqr(constraint_violation);
                            constant_forces.push_back(constraint);
                            RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Penalty_Force_Derivative(s,threshold,std::get<1>(constant_memory),relative_position*direction,offset,spin,force_terms);}
                        right_hand_side.template block<t+d,1>(s*(t+d),0)-=force_direction*right_hand_side_force;}}}}}
    constraint_right_hand_side.resize(rhs.size(),1);
    for(int i=0;i<rhs.size();i++){constraint_right_hand_side(i,0)=rhs[i];}
    system.Flatten_Jacobian_Block(data,force,*this,*rigid_data,terms);
    system.Flatten_Jacobian_Block(data,force,*rigid_data,*this,forces);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(WALL_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(WALL_CONSTRAINT,void)
{
    auto wall_constraint=simulation.force.template Find_Or_Create<WALL_CONSTRAINT<TV>>();
    Parse_Vector(node["lower"],wall_constraint->walls[0]);
    Parse_Vector(node["upper"],wall_constraint->walls[1]);
    return 0;
}
