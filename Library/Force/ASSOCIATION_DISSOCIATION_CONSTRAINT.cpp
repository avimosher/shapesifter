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
        information->value.template block<1+t,1>((1+t)*i,0)=force_memory[information->constraints[i]].second;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<const STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        force_memory[information->constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,information->value.template block<1+t,1>(i*(1+t),0));
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
            memory.second+=increment*information->value.template block<1+t,1>(i*(1+t),0);
        }
        else{
            force_memory[information->constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,increment*information->value.template block<1+t,1>(i*(1+t),0));
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
        for(int i=0;i<interactions.size();i++){
            auto interaction=interactions[i];
            auto structure1=rigid_data->structures[interaction.s1];
            auto structure2=rigid_data->structures[interaction.s2];
            std::pair<int,FORCE_VECTOR> remembered=force_memory[i];
            bool constraint_active=false;
            if(remembered.first==call_count){ // decide whether constraint gets turned off
                T dissociation_rate=1/interaction.base_dissociation_time;
                T cumulative_distribution=1-exp(-dissociation_rate*dt);
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
                T cumulative_distribution=1-exp(-association_rate*dt);
                constraint_active=random.Uniform((T)0,(T)1)<cumulative_distribution;
            }
            if(constraint_active){constraints.push_back(i);}
        }
    }
    //RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
    Matrix<T,t,t+d> angular_to_constraint;angular_to_constraint.setZero();
    angular_to_constraint.template block<t,t>(0,d).setIdentity();
    std::vector<Triplet<LINEAR_CONSTRAINT_MATRIX>> linear_terms;
    std::vector<Triplet<ANGULAR_CONSTRAINT_MATRIX>> angular_terms;
    constraint_rhs.resize(constraints.size()*(1+t));
    stored_forces.resize(constraints.size()*(1+t));
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
        LOG::cout<<"Direction: "<<direction.transpose()<<std::endl;
        auto dC_dX2=RIGID_STRUCTURE_INDEX_MAP<TV>::dC_dA(*structure2,interaction.v2,x1,x2,direction);
        auto dC_dX1=RIGID_STRUCTURE_INDEX_MAP<TV>::dC_dA(*structure1,interaction.v1,x1,x2,direction);
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,interaction.s2,dC_dX2));
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,interaction.s1,-dC_dX1));
        rhs(0,0)=-(direction.norm()-(T).2);


        //auto dC_dX1=index_map.Velocity_Map(*structure1,interaction.v1);
        //auto dC_dX2=index_map.Velocity_Map(*structure2,interaction.v2);
        //linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,interaction.s2,dC_dX2));
        //linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,interaction.s1,-dC_dX1));
        //rhs.template block<d,1>(0,0)=-displacement;

        ROTATION<TV> R1_current=ROTATION<TV>::From_Rotation_Vector(structure1->twist.angular);
        ROTATION<TV> R2_current=ROTATION<TV>::From_Rotation_Vector(structure2->twist.angular);
        ROTATION<TV> R1_base=(R1_current.inverse()*structure1->frame.orientation)*interaction.relative_orientation.inverse();
        ROTATION<TV> R2_base=R2_current.inverse()*structure2->frame.orientation;
        ROTATION<TV> RC=Find_Appropriate_Rotation(R1_current*R1_base,R2_current*R2_base);
        T_SPIN first_rotation_error_vector,second_rotation_error_vector;
        /*for(int j=0;j<t;j++){
            CONSTRAINT_VECTOR constraint1;constraint1.setZero();
            constraint1.template block<1,t>(0,d)=dCdR1.template block<1,t>(j,0);
            CONSTRAINT_VECTOR constraint2;constraint2.setZero();
            constraint2.template block<1,t>(0,d)=dCdR2.template block<1,t>(j,0);
            int index=constraints.size()*d+i*t+j;
            terms.push_back(Triplet<CONSTRAINT_VECTOR>(index,interaction.s1,constraint1));
            terms.push_back(Triplet<CONSTRAINT_VECTOR>(index,interaction.s2,-constraint2));
            }*/
        
        ANGULAR_CONSTRAINT_MATRIX dC_dA1=RIGID_STRUCTURE_INDEX_MAP<TV>::Construct_Constraint_Matrix(R1_current,R1_base*RC,first_rotation_error_vector)*angular_to_constraint;
        ANGULAR_CONSTRAINT_MATRIX dC_dA2=RIGID_STRUCTURE_INDEX_MAP<TV>::Construct_Constraint_Matrix(R2_current,R2_base*RC,second_rotation_error_vector)*angular_to_constraint;
        /*for(int ii=0;ii<dC_dA1.rows();ii++){
            for(int jj=0;jj<dC_dA1.cols();jj++){
                LOG::cout<<"("<<ii<<","<<jj<<"): "<<dC_dA1(ii,jj)<<std::endl;
                }}*/
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,interaction.s2,dC_dA2.eval()));
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,interaction.s1,-dC_dA1.eval()));
        auto total_rotation_error=second_rotation_error_vector-first_rotation_error_vector;
        //LOG::cout<<"Total rotation error: "<<total_rotation_error.transpose()<<std::endl;
        rhs.template block<t,1>(1,0)=-total_rotation_error;
        constraint_rhs.template block<1+t,1>((1+t)*i,0)=rhs;
        auto& remembered=force_memory[interaction_index];
        if(remembered.first!=call_count){remembered.second.setZero();}
        std::cout<<"First memory: "<<remembered.second.transpose()<<std::endl;
        stored_forces.template block<1+t,1>((1+t)*i,0)=remembered.second;
        T right_hand_force=remembered.second[0];
        T_SPIN right_hand_torque=remembered.second.template block<t,1>(1,0);
        /*LOG::cout<<"Force: "<<right_hand_force<<" torque: "<<right_hand_torque.transpose()<<std::endl;
        LOG::cout<<"dC_dA1: "<<dC_dA1<<std::endl;
        LOG::cout<<"dC_dA2: "<<dC_dA2<<std::endl;
        LOG::cout<<"dC_dX1: "<<dC_dX1<<std::endl;
        LOG::cout<<"dC_dX2: "<<dC_dX2<<std::endl;
        LOG::cout<<"RH contribution 1: "<<(dC_dA1.transpose()*right_hand_torque+dC_dX1.transpose()*right_hand_force).transpose()<<std::endl;
        LOG::cout<<"RH contribution 2: "<<-(dC_dA2.transpose()*right_hand_torque+dC_dX2.transpose()*right_hand_force).transpose()<<std::endl;*/
        right_hand_side.template block<t+d,1>(interaction.s1*(t+d),0)+=dC_dA1.transpose()*right_hand_torque+dC_dX1.transpose()*right_hand_force;
        right_hand_side.template block<t+d,1>(interaction.s2*(t+d),0)-=dC_dA2.transpose()*right_hand_torque+dC_dX2.transpose()*right_hand_force;
    }
    constraint_terms.resize(constraints.size()*(1+t),rigid_data->Velocity_DOF());
    Flatten_Matrices(linear_terms,constraints.size(),angular_terms,constraint_terms);
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
