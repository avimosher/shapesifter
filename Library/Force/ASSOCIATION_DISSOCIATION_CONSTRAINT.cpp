#include <Force/ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.find("RIGID_STRUCTURE_DATA")->second);
    for(int i=0;i<interactions.size();i++){
        auto interaction=interactions[i];
        auto rigid_structure1=rigid_data->structures[interaction.s1];
        auto rigid_structure2=rigid_data->structures[interaction.s2];
        std::pair<int,T> remembered=force_memory[i];
        if(stochastic){
            if(remembered.first==call_count){ // decide whether constraint gets turned off
                T dissociation_rate=1/interaction.base_dissociation_time;
                T cumulative_distribution=1-exp(-dissociation_rate*remembered_dt);
                constraint_active=random.Uniform((T)0,(T)1)>cumulative_distribution;
            }
            else{ // decide whether constraint gets turned on
                
            }
        }

    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(ASSOCIATION_DISSOCIATION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(ASSOCIATION_DISSOCIATION_CONSTRAINT,void)
{
    auto constraint=std::make_shared<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>();
    
    brownian_force->temperature=node["temperature"].asDouble();
    simulation.force.push_back(constraint);
    return 0;
}
