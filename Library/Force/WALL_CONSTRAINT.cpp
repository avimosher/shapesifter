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
        // check against each wall
        for(int axis=0;axis<d;axis++){
            if(walls.maximum_corner[axis]){
                T constraint_violation=data.domain.maximum_corner[axis]-maxima[axis]-structure->collision_radius;
                if(constraint_violation<0){ // do something about it
                    CONSTRAINT constraint(s,axis,wall);
                    auto& memory=force_memory[constraint];
                    if(constraint_violation<slack_distance){
                        if(std::get<0>(memory)!=call_count){
                            std::get<1>(memory)=0;
                            if(std::get<0>(constant_memory)==call_count){
                                std::get<1>(memory)=std::get(constant_memory)*sqr(constraint_violation);}}
                        
                    }
                }
            }
            if(walls.minimum_corner[axis]){
                T constraint_violation=minima[axis]-data.domain.minimum_corner[axis]-structure->collision_radius;
                if(constraint_violation<0){ // do something about it

                }
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(WALL_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(WALL_CONSTRAINT,void)
{
    auto wall_constraint=simulation.force.template Find_Or_Create<WALL_CONSTRAINT<TV>>();
    return 0;
}
#endif
