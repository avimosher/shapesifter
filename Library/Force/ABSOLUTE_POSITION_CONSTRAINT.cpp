#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/MATRIX_BUNDLE.h>
#include <Force/FORCE.h>
#include <Force/ABSOLUTE_POSITION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/LOG.h>
#include <Utilities/MATH.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <stdexcept>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class T> ROTATION<Matrix<T,3,1>> Find_Appropriate_Rotation(const ROTATION<Matrix<T,3,1>>& rotation1,const ROTATION<Matrix<T,3,1>>& rotation2)
{
    return rotation1.inverse()*ROTATION<Matrix<T,3,1>>((rotation2*rotation1.inverse()).inverse()).Scale_Angle((T).5);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ABSOLUTE_POSITION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    throw std::logic_error("Implementation out of date");
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    Matrix<T,Dynamic,1>& right_hand_side=system.RHS(data,force,*rigid_data);
    Matrix<T,Dynamic,1>& constraint_right_hand_side=system.RHS(data,force,*this);

    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<Triplet<FORCE_VECTOR>> forces;
    constraint_right_hand_side.resize(DOF());
    for(int i=0;i<linear_constraints.size();i++){
        const LINEAR_CONSTRAINT& constraint=linear_constraints[i];
        auto structure=rigid_data->structures[constraint.s];
        constraint_right_hand_side[i]=constraint.magnitude-constraint.direction.dot(structure->frame.position);
        right_hand_side.template block<d,1>(constraint.s*(t+d),0)-=constraint.direction*stored_forces[i];
        CONSTRAINT_VECTOR constraint_vector;constraint_vector.setZero();
        constraint_vector.template block<1,d>(0,0)=constraint.direction.transpose();
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,constraint.s,constraint_vector));
        forces.push_back(Triplet<FORCE_VECTOR>(constraint.s,i,constraint_vector.transpose()));
        //RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Constraint_Force_Derivative(constraint.s,stored_forces[i],constraint_right_hand_side[i]*constraint.direction,TV(),T_SPIN(),force_terms); // TODO: this constraint has no spin dependence!
    }
    for(int i=0;i<angular_constraints.size();i++){
        const ANGULAR_CONSTRAINT& constraint=angular_constraints[i];
        auto structure=rigid_data->structures[constraint.s];
        FRAME<TV> frame=structure->frame;
        ROTATION<TV> current(ROTATION<TV>::From_Rotation_Vector(structure->twist.angular));
        ROTATION<TV> base(current.inverse()*frame.orientation);
        ROTATION<TV> RC=Find_Appropriate_Rotation(frame.orientation,constraint.orientation);
        T_SPIN rotation_error_vector;
        Matrix<T,t,t> dCdR=RIGID_STRUCTURE_INDEX_MAP<TV>::Relative_Orientation_Constraint_Matrix(current,base*RC,constraint.orientation*RC,rotation_error_vector);
        T_SPIN stored_torque;
        for(int j=0;j<t;j++){
            CONSTRAINT_VECTOR constraint_vector;constraint_vector.setZero();
            constraint_vector.template block<1,t>(0,d)=dCdR.template block<1,t>(j,0);
            int index=linear_constraints.size()+i*t+j;
            terms.push_back(Triplet<CONSTRAINT_VECTOR>(index,constraint.s,constraint_vector));
            forces.push_back(Triplet<FORCE_VECTOR>(constraint.s,index,FORCE_VECTOR::Unit(d+j)));
            constraint_right_hand_side[index]=-rotation_error_vector[j];
            stored_torque[j]=stored_forces[index];}
        right_hand_side.template block<t,1>(constraint.s*(t+d)+d,0)-=stored_torque;}
    stored_forces.resize(DOF());
    system.Flatten_Jacobian_Block(data,force,*this,*rigid_data,terms);
    system.Flatten_Jacobian_Block(data,force,*rigid_data,*this,forces);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(ABSOLUTE_POSITION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(ABSOLUTE_POSITION_CONSTRAINT,void)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto absolute_position_constraint=simulation.force.template Find_Or_Create<ABSOLUTE_POSITION_CONSTRAINT<TV>>();
    Json::Value constraints=node["constraints"];
    for(Json::ValueIterator it=constraints.begin();it!=constraints.end();it++){
        int s=rigid_data->Structure_Index((*it)["structure"].asString());
        std::string constraint_type;Parse_String((*it)["type"],constraint_type);
        if(constraint_type=="linear"){
            typename ABSOLUTE_POSITION_CONSTRAINT<TV>::LINEAR_CONSTRAINT constraint;
            constraint.s=s;
            Parse_Vector((*it)["direction"],constraint.direction);
            Parse_Scalar((*it)["magnitude"],constraint.magnitude);
            absolute_position_constraint->linear_constraints.push_back(constraint);}
        else if(constraint_type=="angular"){
            typename ABSOLUTE_POSITION_CONSTRAINT<TV>::ANGULAR_CONSTRAINT constraint;
            constraint.s=s;
            Parse_Rotation((*it)["orientation"],constraint.orientation);
            absolute_position_constraint->angular_constraints.push_back(constraint);}}
    absolute_position_constraint->stored_forces.resize(absolute_position_constraint->DOF());
    absolute_position_constraint->stored_forces.setZero();
    return 0;
}
