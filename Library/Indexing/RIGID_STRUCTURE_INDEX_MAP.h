#ifndef __RIGID_STRUCTURE_INDEX_MAP__
#define __RIGID_STRUCTURE_INDEX_MAP__

#include <Data/FRAME.h>
#include <Data/ROTATION.h>
#include <Utilities/EIGEN_HELPERS.h>

namespace Mechanics{
template<class TV>
class RIGID_STRUCTURE_INDEX_MAP
{
    typedef typename TV::Scalar T;
    typedef Matrix<T,1,TV::RowsAtCompileTime> TV_T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{STATIC_SIZE=TWIST<TV>::STATIC_SIZE,d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};

    RIGID_STRUCTURE_INDEX_MAP(){}
    ~RIGID_STRUCTURE_INDEX_MAP(){}

    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Map_Twist_To_Velocity(const RIGID_STRUCTURE<TV>& structure,const TV& offset){
        return Map_Twist_To_Velocity(structure.frame.orientation*offset);
    }

    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Map_Twist_To_Velocity(const TV& offset){
        Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> unknown_map;
        unknown_map.template block<d,d>(0,0).setIdentity();
        unknown_map.template block<d,ROTATION<TV>::TwistSize>(0,d)=Cross_Product_Matrix(offset).transpose();
        return unknown_map;
    }

    static Matrix<T,d,t> dRotatedOffset_dSpin(const T_SPIN& spin,const TV& offset){
        static const T eps=1e-8;
        T norm_spin=spin.norm();
        TV_T dw_dspin=-sinc(norm_spin/2)/4*spin.transpose();        
        TV spin_normspin=(norm_spin>eps)?(TV)(spin/norm_spin):TV::UnitX();
        Matrix<T,3,3> dq_dspin=cos(norm_spin/2)/2*spin_normspin*spin_normspin.transpose()+sinc(norm_spin/2)/2*(Matrix<T,t,t>::Identity()-spin_normspin*spin_normspin.transpose());
        T w=cos(norm_spin/2);
        TV q=sinc(norm_spin/2)*spin/2;
        return -2*Cross_Product_Matrix(offset)*(q*dw_dspin+w*dq_dspin);
    }

    static Matrix<T,1,t+d> dXN_dA(const RIGID_STRUCTURE<TV>& structure,const TV& object_offset,const TV& x1,const TV& x2){
        // loosely, (x2-x1)^T/|x2-x1|*(dx2/dx-dx1/dx)
        // should be t+d columns, one row
        static const T eps=1e-8;
        ROTATION<TV> orientation=ROTATION<TV>::From_Rotation_Vector(structure.twist.angular).inverse()*structure.frame.orientation;

        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=dRotatedOffset_dSpin(structure.twist.angular,orientation*object_offset);//-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);
        TV delta=x2-x1;
        T delta_norm=delta.norm();
        TV delta_over_norm=(delta_norm>eps?(TV)(delta/delta_norm):TV::UnitX());
        return delta_over_norm.transpose()*dx_da;
    }

    static Matrix<T,t+d,t+d> dForce_dTwist(const RIGID_STRUCTURE<TV>& structure,const TV& object_offset,const TV& x1,const TV& x2,const TV& direction){
        TV current_offset=structure.frame.orientation*object_offset;
        ROTATION<TV> delta_orientation=ROTATION<TV>::From_Rotation_Vector(structure.twist.angular).inverse();
        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=dRotatedOffset_dSpin(structure.twist.angular,delta_orientation*current_offset);

        T distance=std::max((T)1e-3,(x2-x1).norm());
        Matrix<T,d,t+d> dd_da=dx_da/distance-direction/cube(distance)*direction.transpose()*dx_da;

        TV d_n=direction.normalized();
        Matrix<T,d,t+d> dr_da;
        dr_da.template block<d,d>(0,0).setZero();
        dr_da.template block<d,t>(0,d)=dRotatedOffset_dSpin(structure.twist.angular,delta_orientation*current_offset);

        Matrix<T,t+d,t+d> dF_da;
        dF_da.template block<d,t+d>(0,0)=dd_da;
        dF_da.template block<t,t+d>(d,0)=-Cross_Product_Matrix(d_n)*dr_da+Cross_Product_Matrix(current_offset)*dd_da;
        return dF_da;
    }

    static Matrix<T,1,t+d> dConstraint_dTwist(const RIGID_STRUCTURE<TV>& structure,const TV& object_offset,const TV& x1,const TV& x2,const TV& direction){
        ROTATION<TV> orientation=ROTATION<TV>::From_Rotation_Vector(structure.twist.angular).inverse()*structure.frame.orientation;
        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=dRotatedOffset_dSpin(structure.twist.angular,orientation*object_offset);
        auto distance=std::max((T)1e-3,(x2-x1).norm());
        TV normalized_direction=direction.normalized();
        auto dd_da=dx_da/distance-normalized_direction/distance*normalized_direction.transpose()*dx_da;
        auto final=direction.transpose()*dd_da+normalized_direction.transpose()*(dx_da);
        return final;
    }

    static Matrix<T,d,d> dForce_dVelocity(const TV& relative_position,int s1,int s2){
        T distance=std::max((T)1e-3,relative_position.norm());
        int sign=s1==s2?1:-1;
        return sign*(Matrix<T,d,d>::Identity()/distance-relative_position/cube(distance)*relative_position.transpose());
    }

    static Matrix<T,d,d> dForce_dSpin(const TV& relative_position,int s1,int s2,const T_SPIN& spin,const TV& rotated_offset){
        return dForce_dVelocity(relative_position,s1,s2)*dRotatedOffset_dSpin(spin,rotated_offset);
    }

    static Matrix<T,t,d> dTorque_dVelocity(const TV& relative_position,int s1,int s2,const TV& rotated_offset){
        return Cross_Product_Matrix(rotated_offset)*dForce_dVelocity(relative_position,s1,s2);
    }
    
    static Matrix<T,t,t> dTorque_dSpin(const TV& relative_position,int s1,int s2,const T_SPIN& spin,const TV& rotated_offset){
        T distance=std::max((T)1e-3,relative_position.norm());        
        return Cross_Product_Matrix(dRotatedOffset_dSpin(spin,rotated_offset))*relative_position/distance+Cross_Product_Matrix(rotated_offset)*dForce_dSpin(relative_position,s1,s2,spin,rotated_offset);
    }

    static Matrix<T,3,3> Compute_Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,const int composed_rotation_sign)
    {
        auto orientation=rotation.Rotation_Vector();
        auto angle=orientation.norm();TV axis=(fabs(angle)>1e-8?orientation.normalized():TV::UnitX());
        T s=sin(angle/2);T s_over_angle=sinc(angle/2)/2,c=cos(angle/2);
        ROTATION<TV> composed_rotation=rotation*relative_rotation;
    
        Matrix<T,3,3> axis_projection=axis*axis.transpose();
        Matrix<T,3,3> axis_orthogonal_projection=Matrix<T,3,3>::Identity()-axis_projection;
        TV relative_rotation_vec=relative_rotation.vec();
        Matrix<T,3,3> relative_rotation_cross_product_matrix=Cross_Product_Matrix(relative_rotation_vec);

        Matrix<T,3,3> dudw=(c/2)*axis_projection+s_over_angle*axis_orthogonal_projection;
        TV dadw=-s/2*axis;
        Matrix<T,3,3> dCdu=Matrix<T,3,3>::Identity()*relative_rotation.w()-relative_rotation_cross_product_matrix;
        TV dCda=relative_rotation.vec();
        return composed_rotation_sign*(dCda*dadw.transpose()+dCdu*dudw);
    }

    static Matrix<T,3,3> Relative_Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,const ROTATION<TV>& target,TV& rotation_error_vector)
    {
        ROTATION<TV> composed_rotation=rotation*relative_rotation;
        rotation_error_vector=composed_rotation.vec()*composed_rotation.Sign()-target.Sign()*target.vec();
        return Compute_Orientation_Constraint_Matrix(rotation,relative_rotation,composed_rotation.Sign());
    }

    static Matrix<T,3,3> Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,TV& rotation_error_vector)
    {
        ROTATION<TV> composed_rotation=rotation*relative_rotation;
        rotation_error_vector=composed_rotation.vec()*composed_rotation.Sign();
        return Compute_Orientation_Constraint_Matrix(rotation,relative_rotation,composed_rotation.Sign());
    }
};
}
#endif
