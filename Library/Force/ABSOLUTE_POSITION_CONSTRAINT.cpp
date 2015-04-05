#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/ABSOLUTE_POSITION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/LOG.h>
#include <Utilities/MATH.h>
#include <Utilities/OSG_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
#include <osg/Geometry>
#include <osg/Geode>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class T> ROTATION<Matrix<T,3,1>> Find_Appropriate_Rotation(const ROTATION<Matrix<T,3,1>>& rotation1,const ROTATION<Matrix<T,3,1>>& rotation2)
{
    return rotation1.inverse()*ROTATION<Matrix<T,3,1>>((rotation2*rotation1.inverse()).inverse()).Scale_Angle((T).5);
}
///////////////////////////////////////////////////////////////////////
template<class T> Matrix<T,3,3> Construct_Constraint_Matrix(const ROTATION<Matrix<T,3,1>>& rotation,const ROTATION<Matrix<T,3,1>>& relative_rotation,const ROTATION<Matrix<T,3,1>>& target,Matrix<T,3,1>& rotation_error_vector)
{
    auto orientation=rotation.Rotation_Vector();
    auto axis=orientation.normalized();auto angle=orientation.norm();
    T s=sin(angle/2);T s_over_angle=sinc(angle/2)/2,c=cos(angle/2);
    ROTATION<Matrix<T,3,1>> rotation_error=rotation*relative_rotation;
    T at=sgn(rotation_error.w());
    rotation_error_vector=rotation_error.vec()*at-sgn(target.w())*target.vec();
    
    auto axis_projection=axis*axis.transpose();
    auto axis_orthogonal_projection=Matrix<T,3,3>::Identity()-axis_projection;
    Matrix<T,3,1> relative_rotation_vec=relative_rotation.vec();
    auto relative_rotation_cross_product_matrix=Cross_Product_Matrix(relative_rotation_vec);

    auto dudw=(c/2)*axis_projection+s_over_angle*axis_orthogonal_projection;
    auto dadw=-s/2*axis;
    auto dCdu=at*(Matrix<T,3,3>::Identity()*relative_rotation.w()-relative_rotation_cross_product_matrix);
    auto dCda=relative_rotation.vec()*at;
    return dCda*dadw.transpose()+dCdu*dudw;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ABSOLUTE_POSITION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.Find("RIGID_STRUCTURE_DATA"));
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    constraint_rhs.resize(Size());
    for(int i=0;i<linear_constraints.size();i++){
        const LINEAR_CONSTRAINT& constraint=linear_constraints[i];
        auto structure=rigid_data->structures[constraint.s];
        FRAME<TV> frame=structure->frame;
        constraint_rhs[i]=constraint.magnitude-constraint.direction.dot(frame.position);
        LOG::cout<<"absolute constraint: "<<constraint.direction.transpose()<<" rhs: "<<constraint_rhs[i]<<std::endl;
        right_hand_side.template block<d,1>(constraint.s*(t+d),0)-=constraint.direction*stored_forces[i];
        CONSTRAINT_VECTOR constraint_vector;constraint_vector.setZero();
        constraint_vector.template block<1,d>(0,0)=constraint.direction.transpose();
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,constraint.s,constraint_vector));
    }
    for(int i=0;i<angular_constraints.size();i++){
        const ANGULAR_CONSTRAINT& constraint=angular_constraints[i];
        auto structure=rigid_data->structures[constraint.s];
        FRAME<TV> frame=structure->frame;
        ROTATION<TV> current(ROTATION<TV>::From_Rotation_Vector(structure->twist.angular));
        ROTATION<TV> base(current.inverse()*frame.orientation);
        ROTATION<TV> RC=Find_Appropriate_Rotation(frame.orientation,constraint.orientation);
        T_SPIN rotation_error_vector;
        Matrix<T,t,t> dCdR=Construct_Constraint_Matrix(current,base*RC,constraint.orientation*RC,rotation_error_vector);
        for(int j=0;j<d;j++){
            CONSTRAINT_VECTOR constraint_vector;constraint_vector.setZero();
            constraint_vector.template block<1,d>(0,0)=dCdR.template block<t,1>(j,0);
            int index=linear_constraints.size()+i*t+j;
            terms.push_back(Triplet<CONSTRAINT_VECTOR>(index,constraint.s,constraint_vector));
            constraint_rhs[index]=rotation_error_vector[j];
        }
    }
    constraint_terms.resize(Size(),(t+d)*rigid_data->structures.size());
    Flatten_Matrix(terms,constraint_terms);
    stored_forces.resize(Size());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ABSOLUTE_POSITION_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(ABSOLUTE_POSITION_CONSTRAINT)
GENERIC_TYPE_DEFINITION(ABSOLUTE_POSITION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(ABSOLUTE_POSITION_CONSTRAINT,void)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(simulation.data.Find("RIGID_STRUCTURE_DATA"));
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
            absolute_position_constraint->linear_constraints.push_back(constraint);
        }
        else if(constraint_type=="angular"){
            typename ABSOLUTE_POSITION_CONSTRAINT<TV>::ANGULAR_CONSTRAINT constraint;
            constraint.s=s;
            Parse_Rotation((*it)["orientation"],constraint.orientation);
            absolute_position_constraint->angular_constraints.push_back(constraint);
        }
    }
    absolute_position_constraint->stored_forces.resize(absolute_position_constraint->Size());
    absolute_position_constraint->stored_forces.setZero();
    return 0;
}
