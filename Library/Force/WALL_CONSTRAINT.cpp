#if 0
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> std::shared_ptr<FORCE_REFERENCE<typename TV::Scalar>> WALL_CONSTRAINT<TV>::
Create_Stored_Force() const
{
    return std::static_pointer_cast<FORCE_REFERENCE<T>>(std::make_shared<STORED_VOLUME_EXCLUSION_CONSTRAINT<T>>());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    if(stochastic){
        for(auto& memory : force_memory){std::get<1>(memory.second)=(T)0;std::get<0>(memory.second)=-1;}
        for(auto& constant_memory : constant_force_memory){std::get<1>(constant_memory.second)=(T)0;std::get<0>(constant_memory.second)=-1;}}


    T slack_distance=-.005;
    for(int s=0;s<rigid_data->structures.size();s++){
        auto structure=rigid_data->structures[s];
        T_SPIN spin=rigid_structure->twist.angular;
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

            for(int w=0;w<2;w++){
                if(walls[w][axis]){
                    T relative_position;
                    if(w==0){
                        relative_position=axis_extrema[w][axis]+structure->frame.position-data.domain[w][axis];}
                    else{
                        relative_position=data.domain[w][axis]-axis_extrema[w][axis]-structure->frame.position[axis];}

                    T constraint_violation=relative_position-threshold;
                    TV offset=axis_extrema[w];
                    if(constraint_violation<0){ // do something about it
                        CONSTRAINT constraint(s,axis,wall);
                        auto& memory=force_memory[constraint];
                        if(constraint_violation<slack_distance){
                            if(std::get<0>(memory)!=call_count){
                                std::get<1>(memory)=0;
                                if(std::get<0>(constant_memory)==call_count){
                                    std::get<1>(memory)=std::get(constant_memory)*sqr(constraint_violation);}}
                            right_hand_side_force=std::get<1>(memory);
                            terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s,RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(spin,offset,relative_position*TV::Unit(axis))));
                            forces.push_back(Triplet<FORCE_VECTOR>(s,constraints.size(),RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(axis_maximum)*TV::Unit(axis)));
                            RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Constraint_Force_Derivative(s,right_hand_side_force,offset,spin,force_terms);
                        }
                        else if(std::get<0>(memory)==call_count || std::get<0>(constant_memory)==call_count){
                            if(std::get<0>(constant_memory)!=call_count){
                                std::get<1>(constant_memory)=std::get<1>(memory)/sqr(constraint_violation);}
                            right_hand_side_force=std::get<1>(constant_memory)*sqr(constraint_violation);
                            constant_forces.push_back(constraint);
                            RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Penalty_Force_Derivative(s,threshold,std::get<1>(constant_memory),relative_position*TV::Unit(axis),offset,spin,force_terms);}}}}}
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(WALL_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(WALL_CONSTRAINT,void)
{
    auto wall_constraint=simulation.force.template Find_Or_Create<WALL_CONSTRAINT<TV>>();
    return 0;
}
#endif
