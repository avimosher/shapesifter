#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/MATRIX_BUNDLE.h>
#include <Force/FORCE.h>
#include <Force/WALL_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/HASHING.h>
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
    information->constraints=constraint_indices;
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
        std::get<MEMORY_COUNT>(constant_force_memory[constant_force_indices[i]])=call_count;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void WALL_CONSTRAINT<TV>::
Identify_Interactions_And_Compute_Errors(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    Matrix<T,Dynamic,1>& right_hand_side=system.RHS(data,force,*rigid_data);
    Matrix<T,Dynamic,1>& constraint_right_hand_side=system.RHS(data,force,*this);
    if(stochastic){
        for(auto& memory : force_memory){std::get<MEMORY_FORCE>(memory.second)=(T)0;std::get<MEMORY_COUNT>(memory.second)=-1;}
        for(auto& constant_memory : constant_force_memory){std::get<MEMORY_FORCE>(constant_memory.second)=(T)0;std::get<MEMORY_COUNT>(constant_memory.second)=-1;}}

    constraints.clear();
    constraint_indices.clear();
    constant_forces.clear();
    constant_force_indices.clear();
    std::vector<T> rhs;
    T slack_distance=-.005;
    T push_out_distance=1e-8;
    for(int s=0;s<rigid_data->structures.size();s++){
        auto structure=rigid_data->structures[s];
        for(int ss=0;ss<structure->Substructures().size();ss++){
            auto& substructure=structure->Substructure(ss);
            T_SPIN spin=structure->twist.angular;
            std::array<TV,2> extrema={structure->frame.orientation*(-substructure.capsule_extent),structure->frame.orientation*substructure.capsule_extent};
            std::array<TV,2> axis_extrema;
            T threshold=substructure.radius;
            TV position=structure->frame*substructure.offset;
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
                        T relative_position=sgn*(axis_extrema[w][axis]+position[axis]-data.domain[w][axis]);
                        T constraint_violation=relative_position-threshold;
                        TV offset=axis_extrema[w];
                        TV direction=sgn*TV::Unit(axis);
                        if(constraint_violation<0){ // do something about it
                            std::array<int,4> indices={s,ss,axis,w};
                            auto& memory=force_memory[indices];
                            auto& constant_memory=constant_force_memory[indices];
                            T right_hand_side_force=0;
                            // TODO: use new/better/correct use of offset calculators here
                            FORCE_VECTOR force_direction=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(offset).transpose()*direction;
                            if(constraint_violation<slack_distance){ // deep into wall - enforce constraint
                                if(std::get<MEMORY_COUNT>(memory)!=call_count){
                                    std::get<MEMORY_FORCE>(memory)=0;
                                    if(std::get<MEMORY_COUNT>(constant_memory)==call_count){
                                        std::get<MEMORY_FORCE>(memory)=std::get<MEMORY_FORCE>(constant_memory)*sqr(constraint_violation);}}
                                right_hand_side_force=std::get<MEMORY_FORCE>(memory);
                                rhs.push_back(constraint_violation-slack_distance-push_out_distance);
                                constraint_indices.push_back(indices);
                                constraints.push_back(CONSTRAINT(spin,right_hand_side_force,offset,force_direction,relative_position*direction));}
                            else if(std::get<MEMORY_COUNT>(memory)==call_count || std::get<MEMORY_COUNT>(constant_memory)==call_count){
                                if(std::get<MEMORY_COUNT>(constant_memory)!=call_count){
                                    std::get<MEMORY_FORCE>(constant_memory)=std::get<MEMORY_FORCE>(memory)/sqr(constraint_violation);}
                                right_hand_side_force=std::get<MEMORY_FORCE>(constant_memory)*sqr(constraint_violation);
                                constant_force_indices.push_back(indices);
                                constant_forces.push_back(CONSTANT_FORCE(spin,offset,relative_position*direction,threshold));}
                            right_hand_side.template block<t+d,1>(s*(t+d),0)+=force_direction*right_hand_side_force;}}}}}}
    constraint_right_hand_side=Map<Matrix<T,Dynamic,1>>(rhs.data(),rhs.size());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void WALL_CONSTRAINT<TV>::
Compute_Derivatives(DATA<TV>& data,FORCE<TV>& force,MATRIX_BUNDLE<TV>& system)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    std::vector<Triplet<T>>& force_terms=system.Matrix_Block_Terms(data,force,*rigid_data);
    std::vector<Triplet<CONSTRAINT_VECTOR>> constraint_velocity_terms;
    std::vector<Triplet<FORCE_VECTOR>> balance_force_terms;

    for(int i=0;i<constraints.size();i++){
        auto& constraint=constraints[i];
        auto& indices=constraint_indices[i];
        auto& memory=force_memory[indices];
        constraint_velocity_terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,indices[0],RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(constraint.spin,constraint.offset,constraint.relative_position)));
        balance_force_terms.push_back(Triplet<FORCE_VECTOR>(indices[0],i,constraint.force_direction));
        RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Constraint_Force_Derivative(indices[0],constraint.force,constraint.relative_position,constraint.offset,constraint.spin,force_terms);}

    for(int i=0;i<constant_forces.size();i++){
        auto& constant_force=constant_forces[i];
        auto& indices=constant_force_indices[i];
        auto& constant_memory=constant_force_memory[indices];
        RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Penalty_Force_Derivative(indices[0],constant_force.threshold,std::get<MEMORY_FORCE>(constant_memory),constant_force.relative_position,constant_force.offset,constant_force.spin,force_terms);}

    system.Flatten_Jacobian_Block(data,force,*this,*rigid_data,constraint_velocity_terms);
    system.Flatten_Jacobian_Block(data,force,*rigid_data,*this,balance_force_terms);
}
///////////////////////////////////////////////////////////////////////
#ifdef VIEWER
template<class TV> void WALL_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
    osg::Group* group=node->asGroup();
    osg::Group* wall_constraint_group=(osg::Group*)getNamedChild(group,Static_Name());
    if(!wall_constraint_group){
        wall_constraint_group=new osg::Group();
        wall_constraint_group->setName(Static_Name());

        const auto& domain=data.domain;
        osg::Vec4 color(1.0f,1.0f,1.0f,1.0f);
        TV edge_lengths=domain.Edge_Lengths();
        for(int axis=0;axis<d;axis++){
            for(int w=0,sgn=1;w<2;w++,sgn-=2){
                if(walls[w][axis]){
                    int segments=6;
                    for(int i=0;i<=segments;i++){
                        // walk across each dimensions in fractions
                        osg::Geode* horizontal=createLine(color);
                        std::array<TV,2> horizontal_points;
                        for(int j=0;j<2;j++){
                            horizontal_points[j][axis]=domain[w][axis];
                            horizontal_points[j][(axis+1)%d]=domain[0][(axis+1)%d]+i/(T)segments*edge_lengths[(axis+1)%d];
                            horizontal_points[j][(axis+2)%d]=domain[j][(axis+2)%d];}
                        updateLine(horizontal,horizontal_points);
                        osg::Geode* vertical=createLine(color);
                        std::array<TV,2> vertical_points;
                        for(int j=0;j<2;j++){
                            vertical_points[j][axis]=domain[w][axis];
                            vertical_points[j][(axis+2)%d]=domain[0][(axis+2)%d]+i/(T)segments*edge_lengths[(axis+2)%d];
                            vertical_points[j][(axis+1)%d]=domain[j][(axis+1)%d];}
                        updateLine(vertical,vertical_points);
                        wall_constraint_group->addChild(horizontal);
                        wall_constraint_group->addChild(vertical);
                    }}}}
        group->addChild(wall_constraint_group);
    }
}
#endif
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(WALL_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(WALL_CONSTRAINT,void)
{
    auto wall_constraint=simulation.force.template Find_Or_Create<WALL_CONSTRAINT<TV>>();
    Parse_Vector(node["lower"],wall_constraint->walls[0]);
    Parse_Vector(node["upper"],wall_constraint->walls[1]);
    return 0;
}
