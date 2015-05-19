#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Force/ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
#include <Force/FORCE.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/RANDOM.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> std::shared_ptr<FORCE_REFERENCE<typename TV::Scalar>> ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Create_Stored_Force() const
{
    return std::static_pointer_cast<FORCE_REFERENCE<T>>(std::make_shared<STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    information->constraints=constraints;
    information->value.resize(information->Size());
    information->value.setZero();
    for(int i=0;i<information->constraints.size();i++){
        information->value.template block<d+t,1>((d+t)*i,0)=force_memory[information->constraints[i]].second;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<const STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        force_memory[information->constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,information->value.template block<d+t,1>(i*(d+t),0));
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment)
{
    call_count+=increment;
    auto information=std::static_pointer_cast<STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        if(force_memory.find(information->constraints[i])!=force_memory.end()){
            auto& memory=force_memory[information->constraints[i]];
            memory.first=call_count;
            memory.second+=increment*information->value.template block<d+t,1>(i*(d+t),0);
        }
        else{
            force_memory[information->constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,increment*information->value.template block<d+t,1>(i*(d+t),0));
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
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    RANDOM<T> random;
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    if(stochastic){
        for(auto memory : force_memory){memory.second.second.setZero();}
        constraints.clear();

        for(int i=0;i<interaction_types.size();i++){
            // build acceleration structure out of all binding sites
            KdBVH<T,3,int> hierarchy_second_binder;
            
            for(auto first_binder : interaction_types(i).first_binders){
                auto structure1=rigid_data->structures[first_binder.index];
                PROXIMITY_SEARCH<TV> proximity_search(first_binder);
                ROTATION<TV> binder1_frame=structure1->frame*ROTATION<TV>::From_Rotated_Vector(TV::Unit(1),first_binder.second);
                BVIntersect(hierarchy_second_binder,proximity_search);
                for(auto candidate : proximity_search.candidates){
                    auto structure2=rigid_data->structure[candidate.index];
                    TV direction=data.Minimum_Offset(structure1->frame*interaction_type.v1,structure2->frame*interaction_type.v2);
                    T bond_distance=direction.norm();
                    T orientation_compatibility=T(),position_compatibility=T();
                    if(bond_distance<interaction_type.bond_distance_threshold){
                        position_compatibility=1-bond_distance/interaction_type.bond_distance_threshold;
                        ROTATION<TV> composed_rotation(structure2->frame.orientation.inverse()*structure1->frame.orientation*interaction_type.relative_orientation.inverse());
                        orientation_compatibility=std::max((T)0,1-std::abs(composed_rotation.Angle())/interaction_type.bond_orientation_threshold);}
                    T compatibility=orientation_compatibility*position_compatibility;
                    T association_rate=compatibility/interaction_type.base_association_time;
                    T cumulative_distribution=1-exp(-association_rate*dt);
                    constraint_active=random.Uniform((T)0,(T)1)<cumulative_distribution;}
                if(constraint_active){constraints.push_back(SPECIFIED_CONSTRAINT);}
            }
        }
        // acceleration.  really, need to do this for every interaction pair.
        std::vector<AlignedBox<T,d>> bounding_list(rigid_data->Size());
        std::vector<int> index_list(rigid_data->Size());
        for(int s=0;s<rigid_data->Size();s++){
            auto structure=rigid_data->structures[s];
            bounding_list[s]=bounding_box(structure);
            index_list[s]=s;
        }


        KdBVH<T,3,int> hierarchy(index_list.begin(),index_list.end(),bounding_list.begin(),bounding_list.end());

        for(int i=0;i<interactions.size();i++){
            auto interaction=interactions[i];
            auto structure1=rigid_data->structures[interaction.s1];
            auto structure2=rigid_data->structures[interaction.s2];
            std::pair<int,FORCE_VECTOR> remembered=force_memory[i];
            bool constraint_active=false;
            if(remembered.first==call_count){ // decide whether constraint gets turned off
                T dissociation_rate=1/interaction.base_dissociation_time;
                T cumulative_distribution=1-exp(-dissociation_rate*dt);
                constraint_active=random.Uniform((T)0,(T)1)>cumulative_distribution;}
            else{ // decide whether constraint gets turned on
                TV direction=data.Minimum_Offset(structure1->frame*interaction.v1,structure2->frame*interaction.v2);
                T bond_distance=direction.norm();
                T orientation_compatibility=T(),position_compatibility=(T)0;
                if(bond_distance<interaction.bond_distance_threshold){
                    position_compatibility=1-bond_distance/interaction.bond_distance_threshold;
                    ROTATION<TV> composed_rotation(structure2->frame.orientation.inverse()*structure1->frame.orientation*interaction.relative_orientation.inverse());
                    orientation_compatibility=std::max((T)0,1-std::abs(composed_rotation.Angle())/interaction.bond_orientation_threshold);}
                T compatibility=orientation_compatibility*position_compatibility;
                T association_rate=compatibility/interaction.base_association_time;
                T cumulative_distribution=1-exp(-association_rate*dt);
                constraint_active=random.Uniform((T)0,(T)1)<cumulative_distribution;}
            if(constraint_active){constraints.push_back(i);}}}

    Matrix<T,t,t+d> angular_to_constraint;angular_to_constraint.setZero();
    angular_to_constraint.template block<t,t>(0,d).setIdentity();
    std::vector<Triplet<LINEAR_CONSTRAINT_MATRIX>> linear_terms;
    std::vector<Triplet<ANGULAR_CONSTRAINT_MATRIX>> angular_terms;
    constraint_rhs.resize(constraints.size()*(d+t));
    stored_forces.resize(constraints.size()*(d+t));
    for(int i=0;i<constraints.size();i++){
        int interaction_index=constraints[i];
        auto interaction=interactions[interaction_index];
        auto structure1=rigid_data->structures[interaction.s1];
        auto structure2=rigid_data->structures[interaction.s2];
        FORCE_VECTOR rhs;

        auto x1=structure1->frame*interaction.v1;
        auto x2=structure2->frame*interaction.v2;
        //TV direction=structure1->Displacement(data,*structure2,offset1,offset2).normalized(); // use core-core direction for stability reasons
        TV direction=data.Minimum_Offset(x1,x2); // can't easily use point-point distance then, though.
        
        LINEAR_CONSTRAINT_MATRIX dC_dX1=RIGID_STRUCTURE_INDEX_MAP<TV>::Velocity_Map(*structure1,interaction.v1);
        LINEAR_CONSTRAINT_MATRIX dC_dX2=RIGID_STRUCTURE_INDEX_MAP<TV>::Velocity_Map(*structure2,interaction.v2);
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,interaction.s2,dC_dX2));
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,interaction.s1,-dC_dX1));
        rhs.template block<d,1>(0,0)=-direction;

        ROTATION<TV> R1_current=ROTATION<TV>::From_Rotation_Vector(structure1->twist.angular);
        ROTATION<TV> R2_current=ROTATION<TV>::From_Rotation_Vector(structure2->twist.angular);
        ROTATION<TV> R1_base=(R1_current.inverse()*structure1->frame.orientation)*interaction.relative_orientation.inverse();
        ROTATION<TV> R2_base=R2_current.inverse()*structure2->frame.orientation;
        ROTATION<TV> RC=Find_Appropriate_Rotation(R1_current*R1_base,R2_current*R2_base);
        T_SPIN first_rotation_error_vector,second_rotation_error_vector;
        
        ANGULAR_CONSTRAINT_MATRIX dC_dA1=RIGID_STRUCTURE_INDEX_MAP<TV>::Construct_Constraint_Matrix(R1_current,R1_base*RC,first_rotation_error_vector)*angular_to_constraint;
        ANGULAR_CONSTRAINT_MATRIX dC_dA2=RIGID_STRUCTURE_INDEX_MAP<TV>::Construct_Constraint_Matrix(R2_current,R2_base*RC,second_rotation_error_vector)*angular_to_constraint;
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,interaction.s2,dC_dA2.eval()));
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,interaction.s1,-dC_dA1.eval()));
        T_SPIN total_rotation_error=second_rotation_error_vector-first_rotation_error_vector;
        rhs.template block<t,1>(d,0)=-total_rotation_error;
        constraint_rhs.template block<d+t,1>((d+t)*i,0)=rhs;
        auto& remembered=force_memory[interaction_index];
        if(remembered.first!=call_count){remembered.second.setZero();}
        stored_forces.template block<d+t,1>((d+t)*i,0)=remembered.second;
        TV right_hand_force=remembered.second.template block<d,1>(0,0);
        T_SPIN right_hand_torque=remembered.second.template block<t,1>(d,0);
        right_hand_side.template block<d+t,1>(interaction.s1*(d+t),0)+=dC_dA1.transpose()*right_hand_torque+dC_dX1.transpose()*right_hand_force;
        right_hand_side.template block<d+t,1>(interaction.s2*(d+t),0)-=dC_dA2.transpose()*right_hand_torque+dC_dX2.transpose()*right_hand_force;
    }
    constraint_terms.resize(constraints.size()*(d+t),rigid_data->Velocity_DOF());
    Flatten_Matrices(linear_terms,d*constraints.size(),angular_terms,constraint_terms);
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(ASSOCIATION_DISSOCIATION_CONSTRAINT)
GENERIC_TYPE_DEFINITION(ASSOCIATION_DISSOCIATION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(ASSOCIATION_DISSOCIATION_CONSTRAINT,void)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
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
