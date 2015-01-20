#include <Force/ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Find_Appropriate_Rotation(const ROTATION<TV>& rotation1,const ROTATION<TV>& rotation2)
{
    return rotation1.inverse()*(rotation2*rotation1.inverse()).inverse().Scale_Angle((T).5);
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,ROTATION<TV>::TwistSize,ROTATION<TV>::TwistSize> ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Construct_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,TV& rotation_error)
{
    TV axis=rotation.Axis();
    T angle=rotation.Angle();
    T s=sin(angle/2);T s_over_angle=sinc(angle/2),c=cos(angle/2);

    ROTATION<TV> error=rotation*relative_rotation;
    T at=sign(rotation_error.q.s);
    rotation_error=error.q.v*at;

    auto axis_projection=axis*axis.transpose();
    auto axis_orthogonal_projection=ROTATION_MATRIX::Identity()-axis_projection;
    auto relative_rotation_cross_product_matrix=Cross_Product_Matrix(relative_rotation.q.v);

    auto dudw=(c/2)*axis_projection+s_over_angle*axis_orthogonal_projection;
    auto dadw=-s/2*axis_vector;
    auto dCdu=at*(ROTATION_MATRIX::Identity()*relative_rotation.q.s-relative_rotation_cross_product_matrix);
    auto dCda=relative_rotation.q.v*at;
    return dCda*dadw.transpose()+dCdu*dudw;
}
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
        }

        for(int c=0;c<constraints.size();c++){
            linear_terms.push_back(Triplet<LINEAR_MATRIX>(j,body_index1,index_map.Velocity_Map(offset2)));
            linear_terms.push_back(Triplet<LINEAR_MATRIX>(j,body_index2,-index_map.Velocity_Map(offset1)));
            linear_rhs.push_back(-direction);

            ROTATION<TV> RC=Find_Appropriate_Rotation(R1_current*R1_base,R2_current*R2_base);
            angular_terms.push_back(Triplet<ROTATION_MATRIX>(j,body_index1,Construct_Constraint_Matrix(R1_current,R1_base*RC,first_rotation_error_vector)));
            angular_terms.push_back(Triplet<ROTATION_MATRIX>(j,body_index2,-Construct_Constraint_Matrix(R2_current,R2_base*RC,second_rotation_error_vector)));
            auto total_rotation_error=first_rotation_error_vector-second_rotation_error_vector;
            angular_rhs.push_back(-total_rotation_error);
            
        }
    }
    Flatten_Matrix(linear_terms,linear_constraint_matrix);
    Flatten_Matrix(angular_terms,angular_constraint_matrix);
    constraint_terms<<linear_constraint_matrix,Zeros(),
        Zeros(),angular_constraint_matrix;
    // TODO: build full constraint matrix out of linear and angular parts
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
