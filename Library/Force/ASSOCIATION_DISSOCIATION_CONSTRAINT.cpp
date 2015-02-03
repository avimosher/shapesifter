#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Force/ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
#include <Force/FORCE.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/RANDOM.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Unpack_Forces(const Matrix<T,Dynamic,1>& forces)
{
    call_count++;
    for(int i=0;i<constraints.size();i++){
        force_memory[constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,forces.template block<FullSize,1>(i*FullSize,0));
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Increment_Forces(const Matrix<T,Dynamic,1>& forces)
{
    call_count++;
    for(int i=0;i<constraints.size();i++){
        if(force_memory.count(constraints[i])){
            auto& memory=force_memory[constraints[i]];
            memory.first=call_count;
            memory.second+=forces.template block<FullSize,1>(i*FullSize,0);
        }
        else{
            force_memory[constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,forces.template block<FullSize,1>(i*FullSize,0));
        }
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> ROTATION<TV> ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Find_Appropriate_Rotation(const ROTATION<TV>& rotation1,const ROTATION<TV>& rotation2)
{
    return ROTATION<TV>(rotation1.inverse()*(rotation2*rotation1.inverse()).inverse()).Scale_Angle((T).5);
}
///////////////////////////////////////////////////////////////////////
template<class T> Matrix<T,3,3> Cross_Product_Adapter(const Matrix<T,3,1>& vec)
{
    return Cross_Product_Matrix(vec);
}
template<class T> Matrix<T,1,1> Cross_Product_Adapter(const Matrix<T,1,1>& vec)
{
    return Matrix<T,1,1>();
}
template<class T> Matrix<T,0,0> Cross_Product_Adapter(const Matrix<T,0,1>& vec)
{
    return Matrix<T,0,0>();
}
template<class TV> Matrix<typename TV::Scalar,ROTATION<TV>::TwistSize,ROTATION<TV>::TwistSize> ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Construct_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,T_SPIN& rotation_error)
{
    auto axis=rotation.Axis();
    T angle=rotation.Angle();
    T s=sin(angle/2);T s_over_angle=sinc(angle/2),c=cos(angle/2);

    ROTATION<TV> error=rotation*relative_rotation;
    T at=error.Sign();
    rotation_error=error.Vec()*at;

    auto axis_projection=axis*axis.transpose();
    auto axis_orthogonal_projection=ROTATION_MATRIX::Identity()-axis_projection;
    auto relative_rotation_cross_product_matrix=Cross_Product_Adapter(relative_rotation.Vec());

    auto dudw=(c/2)*axis_projection+s_over_angle*axis_orthogonal_projection;
    auto dadw=-s/2*axis;
    auto dCdu=at*(ROTATION_MATRIX::Identity()*relative_rotation.W()-relative_rotation_cross_product_matrix);
    auto dCda=relative_rotation.Vec()*at;
    return dCda*dadw.transpose()+dCdu*dudw;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    RANDOM<T> random;
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.Find("RIGID_STRUCTURE_DATA"));
    if(stochastic){
        constraints.clear();
        for(int i=0;i<interactions.size();i++){
            auto interaction=interactions[i];
            auto structure1=rigid_data->structures[interaction.s1];
            auto structure2=rigid_data->structures[interaction.s2];
            std::pair<int,FORCE_VECTOR> remembered=force_memory[i];
            bool constraint_active=false;
            if(remembered.first==call_count){ // decide whether constraint gets turned off
                T dissociation_rate=1/interaction.base_dissociation_time;
                T cumulative_distribution=1-exp(-dissociation_rate*remembered_dt);
                constraint_active=random.Uniform((T)0,(T)1)>cumulative_distribution;
            }
            else{ // decide whether constraint gets turned on
                TV direction=data.Minimum_Offset(structure1->frame*interaction.v1,structure2->frame*interaction.v2);
                T bond_distance=direction.norm();
                T orientation_compatibility=T(),position_compatibility=(T)0;
                if(bond_distance<interaction.bond_distance_threshold){
                    position_compatibility=1-bond_distance/interaction.bond_distance_threshold;
                    ROTATION<TV> composed_rotation(structure2->frame.orientation.inverse()*structure1->frame.orientation*interaction.relative_orientation.inverse());
                    orientation_compatibility=std::max((T)0,1-abs(composed_rotation.Angle())/interaction.bond_orientation_threshold);
                }
                T compatibility=orientation_compatibility*position_compatibility;
                T association_rate=compatibility/interaction.base_association_time;
                T cumulative_distribution=1-exp(-association_rate*remembered_dt);
                constraint_active=random.Uniform((T)0,(T)1)<cumulative_distribution;
            }
            if(constraint_active){constraints.push_back(i);}
        }
    }
    RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
    Matrix<T,LinearSize,FullSize> linear_to_constraint;linear_to_constraint.setZero();
    linear_to_constraint.template block<LinearSize,LinearSize>(0,0).setIdentity();
    Matrix<T,TwistSize,FullSize> angular_to_constraint;angular_to_constraint.setZero();
    angular_to_constraint.template block<TwistSize,TwistSize>(0,LinearSize).setIdentity();
    std::vector<Triplet<LINEAR_CONSTRAINT_MATRIX>> linear_terms;
    std::vector<Triplet<ANGULAR_CONSTRAINT_MATRIX>> angular_terms;
    constraint_rhs.resize(constraints.size()*FullSize);
    stored_forces.resize(constraints.size()*FullSize);
    for(int i=0;i<constraints.size();i++){
        auto interaction=interactions[constraints[i]];
        auto structure1=rigid_data->structures[interaction.s1];
        auto structure2=rigid_data->structures[interaction.s2];
        FORCE_VECTOR rhs;

        TV displacement=data.Minimum_Offset(structure1->frame*interaction.v1,structure2->frame*interaction.v2);
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,interaction.s2,index_map.Velocity_Map(*structure2,interaction.v2)));
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,interaction.s1,-index_map.Velocity_Map(*structure1,interaction.v1)));
        rhs.template block<LinearSize,1>(0,0)=-displacement;

        ROTATION<TV> R1_current=ROTATION<TV>::From_Rotation_Vector(structure1->twist.angular);
        ROTATION<TV> R2_current=ROTATION<TV>::From_Rotation_Vector(structure2->twist.angular);
        ROTATION<TV> R1_base=(R1_current.inverse()*structure1->frame.orientation)*interaction.relative_orientation.inverse();
        ROTATION<TV> R2_base=R2_current.inverse()*structure2->frame.orientation;
        ROTATION<TV> RC=Find_Appropriate_Rotation(R1_current*R1_base,R2_current*R2_base);
        T_SPIN first_rotation_error_vector,second_rotation_error_vector;
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,interaction.s1,Construct_Constraint_Matrix(R1_current,R1_base*RC,first_rotation_error_vector)*angular_to_constraint));
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,interaction.s2,-Construct_Constraint_Matrix(R2_current,R2_base*RC,second_rotation_error_vector)*angular_to_constraint));
        auto total_rotation_error=first_rotation_error_vector-second_rotation_error_vector;
        rhs.template block<TwistSize,1>(LinearSize,0)=-total_rotation_error;
        constraint_rhs.template block<FullSize,1>(FullSize*i,0)=rhs;
        auto remembered=force_memory[i];
        if(remembered.first!=call_count){remembered.second.setZero();}
        stored_forces.template block<FullSize,1>(FullSize*i,0)=remembered.second;
    }
    constraint_terms.resize(constraints.size()*FullSize,rigid_data->structures.size()*FullSize);
    Flatten_Matrices(linear_terms,constraints.size()*LinearSize,angular_terms,constraint_terms);
    remembered_dt=dt;
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(ASSOCIATION_DISSOCIATION_CONSTRAINT)
GENERIC_TYPE_DEFINITION(ASSOCIATION_DISSOCIATION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(ASSOCIATION_DISSOCIATION_CONSTRAINT,void)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(simulation.data.Find("RIGID_STRUCTURE_DATA"));
    auto constraint=std::make_shared<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>();
    Json::Value interactions=node["interactions"];
    for(Json::ValueIterator it=interactions.begin();it!=interactions.end();it++){
        typename ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::INTERACTION interaction;
        interaction.s1=rigid_data->Structure_Index((*it)["structure1"].asString());
        Parse_Vector((*it)["offset1"],interaction.v1);
        interaction.s2=rigid_data->Structure_Index((*it)["structure2"].asString());
        Parse_Vector((*it)["offset2"],interaction.v2);
        interaction.bond_distance_threshold=(*it)["bond_distance_threshold"].asDouble();
        interaction.bond_orientation_threshold=(*it)["bond_orientation_threshold"].asDouble();
        interaction.base_association_time=(*it)["base_association_time"].asDouble();
        interaction.base_dissociation_time=(*it)["base_dissociation_time"].asDouble();
        constraint->interactions.push_back(interaction);
    }
    simulation.force.push_back(constraint);
    return 0;
}
