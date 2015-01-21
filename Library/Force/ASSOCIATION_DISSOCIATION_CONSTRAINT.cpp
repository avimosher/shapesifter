#include <Force/ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.find("RIGID_STRUCTURE_DATA")->second);
    if(stochastic){
        constraints.clear();
        for(int i=0;i<interactions.size();i++){
            auto interaction=interactions[i];
            auto rigid_structure1=rigid_data->structures[interaction.s1];
            auto rigid_structure2=rigid_data->structures[interaction.s2];
            std::pair<int,T> remembered=force_memory[i];
            if(remembered.first==call_count){ // decide whether constraint gets turned off
                T dissociation_rate=1/interaction.base_dissociation_time;
                T cumulative_distribution=1-exp(-dissociation_rate*remembered_dt);
                constraint_active=random.Uniform((T)0,(T)1)>cumulative_distribution;
            }
            else{ // decide whether constraint gets turned on
                TV direction=data.Minimum_Offset(bond_position1,bond_position2);
                T bond_distance=direction.norm();
                T orientation_compatibility=T(),position_compatibility=(T)0;
                if(bond_distance<interaction.bond_distance_threshold){
                    position_compatibility=1-bond_distance/interaction.bond_distance_threshold;
                    ROTATION<TV> composed_rotation=rigid_structure2->Rotation().inverse()*rigid_structure1->Rotation()*interaction.relative_orientation.inverse();
                    orientation_compatibility=max((T)0,1-abs(composed_rotation.Angle())/interaction.bond_orientation_threshold);
                }
                T compatibility=orientation_compatibility*position_compatibility;
                T association_rate=compatibility/interaction.base_association_time;
                T cumulative_distribution=1-exp(-association_rate*remembered_dt);
                constraint_active=random.Uniform((T)0,(T)1)<cumulative_distribution;
            }
            if(constraint_active){constraints.push_back(i);}
        }
    }
    for(int i=0;i<constraints.size();i++){
        auto interaction=interactions[constraints[i]];
        auto structure1=rigid_data->structures[interaction.s1];
        auto structure2=rigid_data->structures[interaction.s2];
        std::pair<int,T> remembered=force_memory[i];
        ROTATION<TV> R1_current=ROTATION<TV>::From_Rotation_Vector(structure1->twist.angular);
        ROTATION<TV> R2_current=ROTATION<TV>::From_Rotation_Vector(structure2->twist.angular);
        ROTATION<TV> R1_base=(R1_current.inverse()*structure1->frame.orientation)*interaction.relative_orientation.inverse();
        ROTATION<TV> R2_base=R2_current.inverse()*structure2->frame.orientation;

        ROTATION<TV> RC=Find_Appropriate_Rotation(R1_current*R1_base,R2_current*R2_base);
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index1,Construct_Constraint_Matrix(R1_current,R1_base*RC,first_rotation_error_vector)));
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index1,-Construct_Constraint_Matrix(R2_current,R2_base*RC,second_rotation_error_vector)));
        auto total_rotation_error=first_rotation_error_vector-second_rotation_error_vector;
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
