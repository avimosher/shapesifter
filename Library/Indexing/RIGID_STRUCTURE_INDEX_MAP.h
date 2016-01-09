#ifndef __RIGID_STRUCTURE_INDEX_MAP__
#define __RIGID_STRUCTURE_INDEX_MAP__

#include <Data/FRAME.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/ROTATION.h>
#include <Data/TWIST.h>
#include <Utilities/EIGEN_HELPERS.h>

namespace Mechanics{
template<class TV>
class RIGID_STRUCTURE_INDEX_MAP
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum DEFINITIONS{STATIC_SIZE=TWIST<TV>::STATIC_SIZE,d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,1,d> TV_T;
    typedef Matrix<T,d,d> M_VxV;
    typedef Matrix<T,d,t> M_VxT;
public:

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
        return 2*q.cross(offset)*dw_dspin-2*(Cross_Product_Matrix(offset)*w+Cross_Product_Matrix(q.cross(offset))+Cross_Product_Matrix(q)*Cross_Product_Matrix(offset))*dq_dspin;
    }

    static Matrix<T,t,t> dOffsetCrossForce_dSpin(const T_SPIN& spin,const TV& spun_offset,const TV& force){
        Matrix<T,3,3> dRdS=dRotatedOffset_dSpin(spin,spun_offset);
        Matrix<T,d,d> dFdS;
        for(int i=0;i<t;i++){dFdS.col(i)=dRdS.col(i).cross(force);}
        return dFdS;
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

    static Matrix<T,t,t> dTorque_dSpin(const TV& relative_position,int s1,int s2,const T_SPIN& spin,const TV& rotated_offset1,const TV& rotated_offset2){
        T distance=std::max((T)1e-8,relative_position.norm());
        Matrix<T,t,t> first_term=Cross_Product_Matrix(rotated_offset1)*dForce_dSpin(relative_position,spin,rotated_offset2);
        if(s1==s2){
            return (s1==0?1:-1)*Cross_Product_Matrix(relative_position)*dRotatedOffset_dSpin(spin,rotated_offset2)/distance+first_term;}
        return first_term;
    }

    static Matrix<T,d,d> dPenaltyForce_dVelocity(const TV& relative_position,const T threshold){
        T distance=std::max((T)1e-8,relative_position.norm());
        T one_over_distance=1/distance;
        T threshold_distance=distance-threshold;
        return (Matrix<T,d,d>::Identity()*one_over_distance-relative_position*relative_position.transpose()*cube(one_over_distance))*sqr(threshold_distance)+2*threshold_distance*relative_position*relative_position.transpose()*sqr(one_over_distance);
    }

    static Matrix<T,t,t> dPenaltyTorque_dSpin(const TV& relative_position,int s1,int s2,const T_SPIN& spin,const TV& rotated_offset1,const TV& rotated_offset2,const T threshold){
        T distance=std::max((T)1e-8,relative_position.norm());
        Matrix<T,t,t> first_term=Cross_Product_Matrix(rotated_offset1)*dPenaltyForce_dVelocity(relative_position,threshold)*dRotatedOffset_dSpin(spin,rotated_offset2);
        if(s1==s2){
            return (s1==0?1:-1)*Cross_Product_Matrix(relative_position)*dRotatedOffset_dSpin(spin,rotated_offset2)/distance*sqr(distance-threshold)+first_term;}
        return first_term;
    }

    static void Compute_Constraint_Force_Derivative(const int index1,const int index2,const int s1,const int s2,const T term_force,const TV& relative_position,const TV& rotated_offset1,const TV& rotated_offset2,const T_SPIN& spin,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index1,index2,0,0,dForce_dVelocity(relative_position)*term_force,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,d,t>(index1,index2,0,1,dForce_dSpin(relative_position,spin,rotated_offset2)*term_force,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,t,d>(index1,index2,1,0,dTorque_dVelocity(relative_position,rotated_offset1)*term_force,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,t,t>(index1,index2,1,1,dTorque_dSpin(relative_position,s1,s2,spin,rotated_offset1,rotated_offset2)*term_force,force_terms);
    }

    static void Compute_Constraint_Force_Derivative(const int index,const T scalar_force,const TV& relative_position,const TV& rotated_offset,const T_SPIN& spin,std::vector<Triplet<T>>& force_terms){
        Compute_Constraint_Force_Derivative(index,index,0,0,scalar_force,relative_position,rotated_offset,rotated_offset,spin,force_terms);
    }

    static void Compute_Constraint_Force_Derivatives(const std::array<int,2>& indices,const T scalar_force,const TV& relative_position,const std::array<TV,2>& rotated_offsets,const std::array<T_SPIN,2>& spins,std::vector<Triplet<T>>& force_terms){
        for(int s1=0;s1<2;s1++){ // structure of the force term
            for(int s2=0;s2<2;s2++){ // structure we're taking derivative with respect to
                Compute_Constraint_Force_Derivative(indices[s1],indices[s2],s1,s2,(s1==s2?1:-1)*scalar_force,relative_position,rotated_offsets[s1],rotated_offsets[s2],spins[s2],force_terms);}}
    }

    static void Compute_Penalty_Force_Derivative(const int index1,const int index2,const int s1,const int s2,const T threshold,const T term_force,const TV& relative_position,const TV& rotated_offset1,const TV& rotated_offset2,const T_SPIN& spin,std::vector<Triplet<T>>& force_terms){
        Matrix<T,d,d> dF_dV=dPenaltyForce_dVelocity(relative_position,threshold)*term_force;
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index1,index2,0,0,dF_dV,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,d,t>(index1,index2,0,1,dF_dV*dRotatedOffset_dSpin(spin,rotated_offset2),force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,t,d>(index1,index2,1,0,Cross_Product_Matrix(rotated_offset1)*dF_dV,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,t,t>(index1,index2,1,1,dPenaltyTorque_dSpin(relative_position,s1,s2,spin,rotated_offset1,rotated_offset2,threshold)*term_force,force_terms);
    }

    static void Compute_Penalty_Force_Derivative(const int index,const T threshold,const T force_constant,const TV& relative_position,const TV& rotated_offset,const T_SPIN& spin,std::vector<Triplet<T>>& force_terms){
        Compute_Penalty_Force_Derivative(index,index,0,0,threshold,force_constant,relative_position,rotated_offset,rotated_offset,spin,force_terms);
    }

    static void Compute_Penalty_Force_Derivatives(const std::array<int,2>& indices,const T threshold,const T force_constant,const TV& relative_position,const std::array<TV,2>& rotated_offsets,const std::array<T_SPIN,2>& spins,std::vector<Triplet<T>>& force_terms){
        for(int s1=0;s1<2;s1++){ // structure of the force term
            for(int s2=0;s2<2;s2++){ // structure we're taking derivative with respect to
                Compute_Penalty_Force_Derivative(indices[s1],indices[s2],s1,s2,threshold,(s1==s2?1:-1)*force_constant,relative_position,rotated_offsets[s1],rotated_offsets[s2],spins[s2],force_terms);}}
    }


    /*static void Second_Derivative_Angular(){
        (-2*ostar*dq_db).contractColumns().outer(dw_da);
        +2q.cross(o).contractColumns()*d2w_dadb;
        2(ostar*
        
        }*/

    /*static void d2ao_da2(){
        (q.cross(o).cross(a)).sum()*d2q_da2(a);
    }

    static Matrix<TV,3,3> d2q_da2(const TV& a){
        T na=a.norm();
        TV dna_da=a/na;
        a_terms-=(T).25*sin(na/2)*dna_da.outer(dna_da)/na;
        a_terms+=(T)1.5*cos(na/2)/na^2*(eye-dna_da^2);
        a_terms-=sin(na/2)*((3/na^3)*(eye-dna_da^2));
        }*/

    static Matrix<T,3,3> df2mf1_dVelocity(int term_sign){
        return term_sign*Matrix<T,3,3>::Identity();
    }

    // d(|f2-f1})/dv
    static Matrix<T,3,1> dnf_dVelocity(const TV& f,const T& nf,int term_sign){
        return term_sign*f/nf;
    }

    // d(1/|f2-f1|)/dv
    static Matrix<T,3,1> dnfinv_dVelocity(const TV& f,const T& nf,int term_sign){
        return -1/(nf*nf)*dnf_dVelocity(f,nf,term_sign);
    }

    // d2(1/|f2-f1|)/dv1/dv2
    static Matrix<T,3,3> d2nfinv_dVelocity2(const TV& f,const T& nf,int ts1,int ts2){
        TV dnf_dv1=dnf_dVelocity(f,nf,ts1);
        TV dnf_dv2=dnf_dVelocity(f,nf,ts2);
        TV dnfinv_dv2=dnfinv_dVelocity(f,nf,ts2);
        const M_VxV& df_dv2=df2mf1_dVelocity(ts2);
        return 2/cube(nf)*dnf_dv1*dnf_dv2.transpose()-1/sqr(nf)*ts1*(df_dv2/nf+f*dnfinv_dv2.transpose());
    }

    // d2(|f2-f1|)/dv1/dv2
    static Matrix<T,3,3> d2n_dVelocity2(const TV& f,int ts1,int ts2){
        T nf=f.norm();
        return (df2mf1_dVelocity(ts2)/nf+f*dnfinv_dVelocity(f,nf,ts2).transpose())*df2mf1_dVelocity(ts1);//+f.sum()/nf*
    }

    static Matrix<T,3,3> df_dVelocity(const TV& f,int ts){
        T nf=f.norm();
        return df2mf1_dVelocity(ts)/nf+f*dnfinv_dVelocity(f,nf,ts).transpose();
    }

    static TensorFixedSize<T,Sizes<3,3,3>> Outer_Product(const Matrix<T,3,3>& m,const Matrix<T,3,1>& v,const std::vector<int>& indices){
        TensorFixedSize<T,Sizes<3,3,3>> tensor;
        Matrix<int,3,1> index;index<<0,0,0;
        for(index[0]=0;index[0]<3;index[0]++){
            for(index[1]=0;index[1]<3;index[1]++){
                for(index[2]=0;index[2]<3;index[2]++){
                    tensor(index[0],index[1],index[2])=m(index(indices[0]),index(indices[1]))*v(index(indices[2]));}}}
        return tensor;
    }

    static TensorFixedSize<T,Sizes<3,3,3>> d2f_dVelocity2(const TV& f,int ts1,int ts2){
        T nf=f.norm();
        Matrix<T,3,3> df_dv1=df2mf1_dVelocity(ts1);
        Matrix<T,3,3> df_dv2=df2mf1_dVelocity(ts2);
        TV dnfinv_dv1=dnfinv_dVelocity(f,nf,ts1);
        TV dnfinv_dv2=dnfinv_dVelocity(f,nf,ts2);
        Matrix<T,3,3> d2nfinv_dv2=d2nfinv_dVelocity2(f,nf,ts1,ts2);
        TensorFixedSize<T,Sizes<3,3,3>> t1,t2,t3;
        t1=Outer_Product(df_dv1,dnfinv_dv2,{2,0,1});
        t2=Outer_Product(df_dv2,dnfinv_dv1,{2,1,0});
        t3=Outer_Product(d2nfinv_dv2,f,{0,1,2});
        return t1+t2+t3;
    }
    
    static void Compute_Constraint_Second_Derivatives(const Matrix<T,Dynamic,1>& force_balance_error,const std::array<int,2>& indices,int constraint_index,const T constraint_error,const T scalar_force,const TV& relative_position,const SparseMatrix<T>& f_scaling,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms);


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

    static Matrix<T,3,3> Compute_Simple_Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& full_rotation,int cross_sign)
    {
        TV orientation=rotation.Rotation_Vector();
        ROTATION<TV> remaining_rotation;
        if(cross_sign==1){remaining_rotation=rotation.inverse()*full_rotation;}
        else{remaining_rotation=full_rotation*rotation;}
        TV u=rotation.vec(),v=remaining_rotation.vec();
        T a=rotation.w(),b=remaining_rotation.w();
        T angle=orientation.norm();TV axis=(fabs(angle)>1e-8?orientation.normalized():TV::UnitX());
        T s=sin(angle/2);T s_over_angle=sinc(angle/2)/2,c=cos(angle/2);
    
        Matrix<T,3,3> axis_projection=axis*axis.transpose();
        Matrix<T,3,3> axis_orthogonal_projection=Matrix<T,3,3>::Identity()-axis_projection;
        Matrix<T,3,3> dudw=(c/2)*axis_projection+s_over_angle*axis_orthogonal_projection;
        TV dadw=-s/2*axis;
        return v*dadw.transpose()+cross_sign*b*dudw-Cross_Product_Matrix(v)*dudw;
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
