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

    static Matrix<T,d,t> dRotatedOffset_dSpin(const T_SPIN& spin,const TV& spun_offset){
        static const T eps=1e-8;
        TV offset=ROTATION<TV>::From_Rotation_Vector(spin).inverse()*spun_offset;
        T norm_spin=spin.norm();
        TV_T dw_dspin=-sinc(norm_spin/2)/4*spin.transpose();        
        TV spin_normspin=(norm_spin>eps)?(TV)(spin/norm_spin):TV::UnitX();
        Matrix<T,3,3> dq_dspin=cos(norm_spin/2)/2*spin_normspin*spin_normspin.transpose()+sinc(norm_spin/2)/2*(Matrix<T,t,t>::Identity()-spin_normspin*spin_normspin.transpose());
        T w=cos(norm_spin/2);
        TV q=sinc(norm_spin/2)*spin/2;
        //return 2*q.cross(offset)*dw_dspin-2*w*Cross_Product_Matrix(offset)*dq_dspin+4*w*offset*dw_dspin+
        //    4*q.dot(offset)*dq_dspin;
        return 2*q.cross(offset)*dw_dspin-2*(Cross_Product_Matrix(offset)*w+Cross_Product_Matrix(q.cross(offset))+Cross_Product_Matrix(q)*Cross_Product_Matrix(offset))*dq_dspin;
    }

    static Matrix<T,1,t+d> dConstraint_dTwist(const TV& spin,const TV& offset,const TV& relative_position){
        static const T eps=1e-8;
        T relative_position_norm=relative_position.norm();
        TV normalized_relative_position=(relative_position_norm>eps?(TV)(relative_position/relative_position_norm):TV::UnitX());
        Matrix<T,1,t+d> final;
        final.template block<1,d>(0,0)=normalized_relative_position.transpose();
        final.template block<1,t>(0,d)=normalized_relative_position.transpose()*dRotatedOffset_dSpin(spin,offset);
        return final;
    }

    static Matrix<T,d,d> dForce_dVelocity(const TV& relative_position){
        T distance=std::max((T)1e-8,relative_position.norm());
        return Matrix<T,d,d>::Identity()/distance-relative_position/cube(distance)*relative_position.transpose();
    }

    static Matrix<T,d,d> dForce_dSpin(const TV& relative_position,const T_SPIN& spin,const TV& rotated_offset){
        return dForce_dVelocity(relative_position)*dRotatedOffset_dSpin(spin,rotated_offset);
    }

    static Matrix<T,t,d> dTorque_dVelocity(const TV& relative_position,const TV& rotated_offset){
        return Cross_Product_Matrix(rotated_offset)*dForce_dVelocity(relative_position);
    }

    static Matrix<T,t,t> dTorque_dSpin(const TV& relative_position,int s1,int s2,const T_SPIN& spin,const std::vector<TV>& rotated_offsets){
        T distance=std::max((T)1e-8,relative_position.norm());
        Matrix<T,t,t> first_term=Cross_Product_Matrix(rotated_offsets[s1])*dForce_dSpin(relative_position,spin,rotated_offsets[s2]);
        if(s1==s2){
            return (s1==0?1:-1)*Cross_Product_Matrix(relative_position)*dRotatedOffset_dSpin(spin,rotated_offsets[s2])/distance+first_term;}
        return first_term;
    }

    static Matrix<T,3,3> Compute_Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,const int composed_rotation_sign)
    {
        TV orientation=rotation.Rotation_Vector();
        T angle=orientation.norm();TV axis=(fabs(angle)>1e-8?orientation.normalized():TV::UnitX());
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
