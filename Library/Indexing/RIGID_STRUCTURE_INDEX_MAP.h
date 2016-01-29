#ifndef __RIGID_STRUCTURE_INDEX_MAP__
#define __RIGID_STRUCTURE_INDEX_MAP__

#include <Data/FRAME.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/ROTATION.h>
#include <Data/TWIST.h>
#include <Utilities/EIGEN_HELPERS.h>

namespace Mechanics{
namespace Dimension{
    enum LINEARITY{LINEAR,ANGULAR};
}

template<class TV>
class RIGID_STRUCTURE_INDEX_MAP
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum DEFINITIONS{STATIC_SIZE=TWIST<TV>::STATIC_SIZE,d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,1,d> TV_T;
    typedef Matrix<T,d,d> M_VxV;
    typedef Matrix<T,d,t> M_VxT;
    typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;
    typedef Dimension::LINEARITY LINEARITY;
    static constexpr T epsilon=1e-8;
public:

    RIGID_STRUCTURE_INDEX_MAP(){}
    ~RIGID_STRUCTURE_INDEX_MAP(){}

    // Standard derivatives

    // d(a/|a|)/da
    static M_VxV da_na_dA(const TV& a,const T na){
        if(na<epsilon){return M_VxV::Zero();}
        T one_na=1/na;
        TV a_na=a*one_na;
        return one_na*(M_VxV::Identity()-a_na*a_na.transpose());
    }

    // d(1/|a|)/da
    static TV dnainv_dA(const TV& a,const T& na){
        return -a/cube(na);
    }

    // d2(a/|a|)/da2
    static T_TENSOR d2a_na_dA2(const TV& a,const T na){
        M_VxV da_na_da=da_na_dA(a,na);
        T one_na=1/na;
        TV a_na=a*one_na;
        TV dnainv_da=dnainv_dA(a,na);
        
        return Outer_Product(M_VxV::Identity()-a_na*a_na.transpose(),dnainv_da,{2,0,1})-
            Outer_Product(da_na_da,a_na,{2,0,1})*(1/na)-Outer_Product(da_na_da,a_na,{1,0,2})*(1/na);
    }

    static M_VxV d2na_dA2(const TV& a,const T na){        
        return M_VxV::Identity()/na+a*dnainv_dA(a,na).transpose();
    }


    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Map_Twist_To_Velocity(const RIGID_STRUCTURE<TV>& structure,const TV& offset){
        return Map_Twist_To_Velocity(structure.frame.orientation*offset);
    }

    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Map_Twist_To_Velocity(const TV& offset){
        Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> unknown_map;
        unknown_map.template block<d,d>(0,0).setIdentity();
        unknown_map.template block<d,ROTATION<TV>::TwistSize>(0,d)=Cross_Product_Matrix(offset).transpose();
        return unknown_map;
    }

    static T_SPIN dw_dSpin(const T_SPIN& spin,const T norm_spin){
        return -sinc(norm_spin/2)/4*spin.transpose();
    }

    static M_VxV d2w_dSpin2(const T_SPIN& spin,const T norm_spin){
        TV dna_da=spin/norm_spin;
        return -(T).25*cos(norm_spin/2)*dna_da*dna_da.transpose()-(T).5*sin(norm_spin/2)*d2na_dA2(spin,norm_spin);
    }

    static M_VxV dq_dSpin(const T_SPIN& spin_normspin,const T norm_spin){
        return cos(norm_spin/2)/2*spin_normspin*spin_normspin.transpose()+sinc(norm_spin/2)/2*(Matrix<T,t,t>::Identity()-spin_normspin*spin_normspin.transpose());
    }

    static T_TENSOR d2q_dSpin2(const T_SPIN& spin,const T ns){
        T_SPIN s_ns=spin/ns;
        T_SPIN dns_ds=s_ns;
        M_VxV ds_ns_ds=da_na_dA(spin,ns);
        M_VxV d2ns_ds2=ds_ns_ds;
        T_TENSOR d2s_ns_ds2=d2a_na_dA2(spin,ns);
        
        T s=sinc(ns/2);
        T c=cos(ns/2);
        return Outer_Product(dns_ds*dns_ds.transpose(),spin,{0,1,2})*(-(T).125*s)+
            (Outer_Product(ds_ns_ds,dns_ds,{2,0,1})+Outer_Product(ds_ns_ds,dns_ds,{2,1,0})+Outer_Product(d2ns_ds2,s_ns,{0,1,2}))*(T).5*c+
            d2s_ns_ds2*(s*ns/2);
    }

    // d2(s*o)/ds2
    static T_TENSOR d2so_dSpin2(const T_SPIN& spin,const TV& o){
        ROTATION<TV> rotation(ROTATION<TV>::From_Rotation_Vector(spin));
        TV q=rotation.vec();
        T w=rotation.w();
        M_VxV ostar=Cross_Product_Matrix(o);
        T ns=spin.norm();
        M_VxV d2w_ds2=d2w_dSpin2(spin,ns);
        T_TENSOR d2q_ds2=d2q_dSpin2(spin,ns);
        TV s_ns=spin/ns;
        M_VxV dq_ds=dq_dSpin(s_ns,ns);
        TV dw_ds=dw_dSpin(spin,ns);
        return Outer_Product(-2*ostar*dq_ds,dw_ds,{2,1,0})+
            Outer_Product(d2w_ds2,2*q.cross(o),{0,1,2})+
            Outer_Product(ostar*dq_ds,dw_ds*(-2),{2,0,1})+Cross_Product(-2*dq_ds,ostar*dq_ds)+Cross_Product(2*ostar*dq_ds,dq_ds)+
            Cross_Product(-2*(w*ostar+Cross_Product_Matrix(q.cross(o))+Cross_Product_Matrix(q)*ostar),d2q_ds2);
    }
    
    template<int DTYPE1,int DTYPE2>
    static typename std::enable_if<DTYPE1==Dimension::ANGULAR && DTYPE2==Dimension::ANGULAR,T_TENSOR>::type d2f_dV2(const std::array<int,2>& signs,const T_SPIN& spin,const TV& offset)
    {
        if(signs[0]==signs[1]){
            return d2so_dSpin2(spin,offset);}
        T_TENSOR t;t.setZero();
        return t;
    }

    template<int DTYPE1,int DTYPE2>
    static typename std::enable_if<DTYPE1==Dimension::LINEAR && DTYPE2==Dimension::LINEAR,T_TENSOR>::type d2f_dV2(const std::array<int,2>& signs,const T_SPIN& spin,const TV& offset)
    {
        T_TENSOR t;t.setZero();
        return t;
    }

    // d(s*o)/ds
    static Matrix<T,d,t> dRotatedOffset_dSpin(const T_SPIN& spin,const TV& offset){
        T norm_spin=spin.norm();
        TV_T dw_dspin=dw_dSpin(spin,norm_spin);
        TV spin_normspin=(norm_spin>epsilon)?(TV)(spin/norm_spin):TV::UnitX();
        M_VxV dq_dspin=dq_dSpin(spin_normspin,norm_spin);
        T w=cos(norm_spin/2);
        TV q=sinc(norm_spin/2)*spin/2;
        return 2*q.cross(offset)*dw_dspin-2*(Cross_Product_Matrix(offset)*w+Cross_Product_Matrix(q.cross(offset))+Cross_Product_Matrix(q)*Cross_Product_Matrix(offset))*dq_dspin;
    }

    static Matrix<T,t,t> dOffsetCrossForce_dSpin(const T_SPIN& spin,const TV& offset,const TV& force){
        M_VxV dRdS=dRotatedOffset_dSpin(spin,offset);
        Matrix<T,d,d> dFdS;
        for(int i=0;i<t;i++){dFdS.col(i)=dRdS.col(i).cross(force);}
        return dFdS;
    }

    static Matrix<T,1,t+d> dConstraint_dTwist(const TV& spin,const TV& offset,const TV& relative_position){
        T relative_position_norm=relative_position.norm();
        TV normalized_relative_position=(relative_position_norm>epsilon?(TV)(relative_position/relative_position_norm):TV::UnitX());
        Matrix<T,1,t+d> final;
        final.template block<1,d>(0,0)=normalized_relative_position.transpose();
        final.template block<1,t>(0,d)=normalized_relative_position.transpose()*dRotatedOffset_dSpin(spin,offset);
        return final;
    }

    static Matrix<T,d,d> dForce_dVelocity(const TV& relative_position){
        T distance=std::max(epsilon,relative_position.norm());
        return Matrix<T,d,d>::Identity()/distance-relative_position/cube(distance)*relative_position.transpose();
    }

    static Matrix<T,d,d> dForce_dSpin(const TV& relative_position,const T_SPIN& spin,const TV& offset){
        return dForce_dVelocity(relative_position)*dRotatedOffset_dSpin(spin,offset);
    }

    static Matrix<T,t,d> dTorque_dVelocity(const TV& relative_position,const TV& offset){
        return Cross_Product_Matrix(offset)*dForce_dVelocity(relative_position);
    }

    static Matrix<T,t,t> dTorque_dSpin(const TV& relative_position,int s1,int s2,const T_SPIN& spin,const TV& offset1,const TV& offset2){
        T distance=std::max(epsilon,relative_position.norm());
        Matrix<T,t,t> first_term=Cross_Product_Matrix(offset1)*dForce_dSpin(relative_position,spin,offset2);
        if(s1==s2){
            return (s1==0?1:-1)*Cross_Product_Matrix(relative_position)*dRotatedOffset_dSpin(spin,offset2)/distance+first_term;}
        return first_term;
    }

    static Matrix<T,d,d> dPenaltyForce_dVelocity(const TV& relative_position,const T threshold){
        T distance=std::max(epsilon,relative_position.norm());
        T one_over_distance=1/distance;
        T threshold_distance=distance-threshold;
        return (Matrix<T,d,d>::Identity()*one_over_distance-relative_position*relative_position.transpose()*cube(one_over_distance))*sqr(threshold_distance)+2*threshold_distance*relative_position*relative_position.transpose()*sqr(one_over_distance);
    }

    static Matrix<T,t,t> dPenaltyTorque_dSpin(const TV& relative_position,int s1,int s2,const T_SPIN& spin,const TV& offset1,const TV& offset2,const T threshold){
        T distance=std::max(epsilon,relative_position.norm());
        Matrix<T,t,t> first_term=Cross_Product_Matrix(offset1)*dPenaltyForce_dVelocity(relative_position,threshold)*dRotatedOffset_dSpin(spin,offset2);
        if(s1==s2){
            return (s1==0?1:-1)*Cross_Product_Matrix(relative_position)*dRotatedOffset_dSpin(spin,offset2)/distance*sqr(distance-threshold)+first_term;}
        return first_term;
    }

    static void Compute_Constraint_Force_Derivative(const int index1,const int index2,const int s1,const int s2,const T term_force,const TV& relative_position,const TV& offset1,const TV& offset2,const T_SPIN& spin,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index1,index2,0,0,dForce_dVelocity(relative_position)*term_force,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,d,t>(index1,index2,0,1,dForce_dSpin(relative_position,spin,offset2)*term_force,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,t,d>(index1,index2,1,0,dTorque_dVelocity(relative_position,offset1)*term_force,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,t,t>(index1,index2,1,1,dTorque_dSpin(relative_position,s1,s2,spin,offset1,offset2)*term_force,force_terms);
    }

    static void Compute_Constraint_Force_Derivative(const int index,const T scalar_force,const TV& relative_position,const TV& offset,const T_SPIN& spin,std::vector<Triplet<T>>& force_terms){
        Compute_Constraint_Force_Derivative(index,index,0,0,scalar_force,relative_position,offset,offset,spin,force_terms);
    }

    static void Compute_Constraint_Force_Derivatives(const std::array<int,2>& indices,const T scalar_force,const TV& relative_position,const std::array<TV,2>& offsets,const std::array<T_SPIN,2>& spins,std::vector<Triplet<T>>& force_terms){
        for(int s1=0;s1<2;s1++){ // structure of the force term
            for(int s2=0;s2<2;s2++){ // structure we're taking derivative with respect to
                Compute_Constraint_Force_Derivative(indices[s1],indices[s2],s1,s2,(s1==s2?1:-1)*scalar_force,relative_position,offsets[s1],offsets[s2],spins[s2],force_terms);}}
    }

    static void Compute_Penalty_Force_Derivative(const int index1,const int index2,const int s1,const int s2,const T threshold,const T term_force,const TV& relative_position,const TV& offset1,const TV& offset2,const T_SPIN& spin,std::vector<Triplet<T>>& force_terms){
        Matrix<T,d,d> dF_dV=dPenaltyForce_dVelocity(relative_position,threshold)*term_force;
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index1,index2,0,0,dF_dV,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,d,t>(index1,index2,0,1,dF_dV*dRotatedOffset_dSpin(spin,offset2),force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,t,d>(index1,index2,1,0,Cross_Product_Matrix(offset1)*dF_dV,force_terms);
        Flatten_Matrix_Term<T,t+d,t+d,t,t>(index1,index2,1,1,dPenaltyTorque_dSpin(relative_position,s1,s2,spin,offset1,offset2,threshold)*term_force,force_terms);
    }

    static void Compute_Penalty_Force_Derivative(const int index,const T threshold,const T force_constant,const TV& relative_position,const TV& offset,const T_SPIN& spin,std::vector<Triplet<T>>& force_terms){
        Compute_Penalty_Force_Derivative(index,index,0,0,threshold,force_constant,relative_position,offset,offset,spin,force_terms);
    }

    static void Compute_Penalty_Force_Derivatives(const std::array<int,2>& indices,const T threshold,const T force_constant,const TV& relative_position,const std::array<TV,2>& offsets,const std::array<T_SPIN,2>& spins,std::vector<Triplet<T>>& force_terms){
        for(int s1=0;s1<2;s1++){ // structure of the force term
            for(int s2=0;s2<2;s2++){ // structure we're taking derivative with respect to
                Compute_Penalty_Force_Derivative(indices[s1],indices[s2],s1,s2,threshold,(s1==s2?1:-1)*force_constant,relative_position,offsets[s1],offsets[s2],spins[s2],force_terms);}}
    }


    // d(f2-f1)/dv
    template<int DTYPE>
    static typename std::enable_if<DTYPE==Dimension::LINEAR,M_VxV>::type df_dVelocity(int term_sign,const T_SPIN& spin,const TV& offset){
        return M_VxV::Identity()*term_sign;
    }

    template<int DTYPE>
    static typename std::enable_if<DTYPE==Dimension::ANGULAR,M_VxV>::type df_dVelocity(int term_sign,const T_SPIN& spin,const TV& offset){
        return M_VxV::Identity()*term_sign*dRotatedOffset_dSpin(spin,offset);
    }

    // d(|f2-f1})/dv
    template<int DTYPE>
    static TV dnf_dVelocity(const TV& f,const T& nf,int term_sign,const T_SPIN& spin,const TV& offset){
        return df_dVelocity<DTYPE>(term_sign,spin,offset).transpose()*f/nf;
    }

    // d(1/|f2-f1|)/dv
    template<int DTYPE>
    static TV dnfinv_dVelocity(const TV& f,const T& nf,int term_sign,const T_SPIN& spin,const TV& offset){
        return df_dVelocity<DTYPE>(term_sign,spin,offset).transpose()*dnainv_dA(f,nf);
    }

    // d2(1/|f2-f1|)/dv1/dv2
    template<int DTYPE1,int DTYPE2>
    static M_VxV d2nfinv_dVelocity2(const TV& f,const T& nf,const std::array<int,2>& signs,const std::array<T_SPIN,2>& spins,const std::array<TV,2>& offsets){
        TV dnf_dv1=dnf_dVelocity<DTYPE1>(f,nf,signs[0],spins[0],offsets[0]);
        TV dnf_dv2=dnf_dVelocity<DTYPE2>(f,nf,signs[1],spins[1],offsets[1]);
        M_VxV d2nf_dv2=d2n_dVelocity2<DTYPE1,DTYPE2>(f,signs,spins,offsets);
        return 2/cube(nf)*dnf_dv1*dnf_dv2.transpose()-1/sqr(nf)*d2nf_dv2;
    }

    static M_VxV Contract(const T_TENSOR& t,const TV& v,const std::array<int,3>& indices){
        M_VxV result;result.setZero();
        std::array<int,3> index{};
        for(index[0]=0;index[0]<3;index[0]++){
            for(index[1]=0;index[1]<3;index[1]++){
                for(index[2]=0;index[2]<3;index[2]++){
                    result(index[1],index[2])+=t(index[indices[0]],index[indices[1]],index[indices[2]])*v(index[0]);
                }}}
        return result;
    }

    // d2(|f2-f1|)/dv1/dv2
    template<int DTYPE1,int DTYPE2>    
    static M_VxV d2n_dVelocity2(const TV& f,const std::array<int,2>& signs,const std::array<T_SPIN,2>& spins,const std::array<TV,2>& offsets){
        T nf=f.norm();
        TV f_nf=f/nf;
        M_VxV df_dv1=df_dVelocity<DTYPE1>(signs[0],spins[0],offsets[0]);
        M_VxV df_dv2=df_dVelocity<DTYPE2>(signs[1],spins[1],offsets[1]);
        T_TENSOR d2f_dv2=d2f_dV2<DTYPE1,DTYPE2>(signs,spins[0],offsets[0]);
        return 1/nf*df_dv1.transpose()*(M_VxV::Identity()-f_nf*f_nf.transpose())*df_dv2+
            Contract(d2f_dv2,f_nf,{1,2,0});
    }

    // d((f2-f1)/|f2-f1|)/dv
    template<int DTYPE>
    static M_VxV df_nf_dVelocity(const TV& f,int ts,const T_SPIN& spin,const TV& offset){
        T nf=f.norm();
        return df_dVelocity<DTYPE>(ts,spin,offset)/nf+f*dnfinv_dVelocity<DTYPE>(f,nf,ts,spin,offset).transpose();
    }

    // assumption: the columns of m1 should go in index 0, columns of m2 in index 1, and cross product results in index 2
    static T_TENSOR Cross_Product(const M_VxV& m1,const M_VxV& m2){
        T_TENSOR tensor;
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                TV cross=m1.col(i).cross(m2.col(j));
                for(int k=0;k<3;k++){
                    tensor(i,j,k)=cross(k);}}}
        return tensor;
    }

    // cross product between cross product matrix and dimension 2 of a tensor
    static T_TENSOR Cross_Product(const M_VxV& m,const T_TENSOR& t){
        T_TENSOR tensor;
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                TV tvec;
                for(int k=0;k<3;k++){tvec[k]=t(i,j,k);}
                TV cross=m*tvec;
                for(int k=0;k<3;k++){
                    tensor(i,j,k)=cross(k);}}}
        return tensor;
    }

    static T_TENSOR Outer_Product(const M_VxV& m,const TV& v,const std::vector<int>& indices){
        T_TENSOR tensor;
        Matrix<int,3,1> index;index<<0,0,0;
        for(index[0]=0;index[0]<3;index[0]++){
            for(index[1]=0;index[1]<3;index[1]++){
                for(index[2]=0;index[2]<3;index[2]++){
                    tensor(index[0],index[1],index[2])=m(index(indices[0]),index(indices[1]))*v(index(indices[2]));}}}
        return tensor;
    }

    // d2((f2-f1)/|f2-f1|)/dv2
    template<int DTYPE1,int DTYPE2>
    static T_TENSOR d2f_nf_dVelocity2(const TV& f,const std::array<int,2>& signs,const std::array<T_SPIN,2>& spins,const std::array<TV,2>& offsets){
        T nf=f.norm();
        M_VxV df_dv1=df_dVelocity<DTYPE1>(signs[0],spins[0],offsets[0]);//checked
        M_VxV df_dv2=df_dVelocity<DTYPE2>(signs[1],spins[1],offsets[1]);//checked
        TV dnfinv_dv1=dnfinv_dVelocity<DTYPE1>(f,nf,signs[0],spins[0],offsets[0]);//checked
        TV dnfinv_dv2=dnfinv_dVelocity<DTYPE2>(f,nf,signs[1],spins[1],offsets[1]);//checked
        M_VxV d2nfinv_dv2=d2nfinv_dVelocity2<DTYPE1,DTYPE2>(f,nf,signs,spins,offsets);//checked
        T_TENSOR d2f_dv2=d2f_dV2<DTYPE1,DTYPE2>(signs,spins[0],offsets[0]);//basically checked
        return Outer_Product(df_dv1,dnfinv_dv2,{2,0,1})+
            Outer_Product(df_dv2,dnfinv_dv1,{2,1,0})+
            Outer_Product(d2nfinv_dv2,f,{0,1,2})+
            d2f_dv2*(1/nf);
    }

    // d(r x (f2-f1)/|f2-f1|)/da
    template<int DTYPE>
    static M_VxV dtau_dA(const TV& f,int sgn,const T_SPIN& spin,const TV& offset){
        M_VxV dr_da=dRotatedOffset_dSpin(spin,offset);
        M_VxV df_da=df_nf_dVelocity<DTYPE>(f,sgn,spin,offset);
        TV spun_offset=ROTATION<TV>::From_Rotation_Vector(spin)*offset;
        return -Cross_Product_Matrix(f.normalized())*dr_da+Cross_Product_Matrix(spun_offset)*df_da;
    }

    // d2(r x (f2-f1)/|f2-f1|)/dv2
    template<int DTYPE1,int DTYPE2>
    static T_TENSOR d2tau_dV2(const TV& f,const TV& r,const std::array<int,2>& signs,const std::array<T_SPIN,2>& spins,const std::array<TV,2>& offsets){
        M_VxV dr_da=dRotatedOffset_dSpin(spins[0],offsets[0]); // checked
        M_VxV dr_db=dRotatedOffset_dSpin(spins[1],offsets[1]); // checked
        M_VxV df_da=df_nf_dVelocity<DTYPE1>(f,signs[0],spins[0],offsets[0]); // checked
        M_VxV df_db=df_nf_dVelocity<DTYPE2>(f,signs[1],spins[1],offsets[1]); // checked
        TV f_nf=f.normalized();
        T_TENSOR d2r_da2=d2f_dV2<DTYPE1,DTYPE2>(signs,spins[0],offsets[0]); //checked
        T_TENSOR d2f_nf_da2=d2f_nf_dVelocity2<DTYPE1,DTYPE2>(f,signs,spins,offsets); //checked
        return Cross_Product(-Cross_Product_Matrix(f_nf),d2r_da2)+
            Cross_Product(dr_da,df_db)+
            Cross_Product(-df_da,dr_db)+ // TODO: correct?
            Cross_Product(Cross_Product_Matrix(r),d2f_nf_da2);
    }
    
    static void Compute_Constraint_Second_Derivatives(const Matrix<T,Dynamic,1>& force_balance_error,const std::array<int,2>& indices,int constraint_index,const T constraint_error,const T scalar_force,const TV& relative_position,const SparseMatrix<T>& f_scaling,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms);


    static M_VxV Compute_Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,const int composed_rotation_sign)
    {
        TV orientation=rotation.Rotation_Vector();
        T angle=orientation.norm();TV axis=(fabs(angle)>1e-8?orientation.normalized():TV::UnitX());
        T s=sin(angle/2);T s_over_angle=sinc(angle/2)/2,c=cos(angle/2);
        ROTATION<TV> composed_rotation=rotation*relative_rotation;
    
        M_VxV axis_projection=axis*axis.transpose();
        M_VxV axis_orthogonal_projection=M_VxV::Identity()-axis_projection;
        TV relative_rotation_vec=relative_rotation.vec();
        M_VxV relative_rotation_cross_product_matrix=Cross_Product_Matrix(relative_rotation_vec);

        M_VxV dudw=(c/2)*axis_projection+s_over_angle*axis_orthogonal_projection;
        TV dadw=-s/2*axis;
        M_VxV dCdu=M_VxV::Identity()*relative_rotation.w()-relative_rotation_cross_product_matrix;
        TV dCda=relative_rotation.vec();
        return composed_rotation_sign*(dCda*dadw.transpose()+dCdu*dudw);
    }

    static M_VxV Compute_Simple_Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& full_rotation,int cross_sign)
    {
        TV orientation=rotation.Rotation_Vector();
        ROTATION<TV> remaining_rotation;
        if(cross_sign==1){remaining_rotation=rotation.inverse()*full_rotation;}
        else{remaining_rotation=full_rotation*rotation;}
        TV u=rotation.vec(),v=remaining_rotation.vec();
        T a=rotation.w(),b=remaining_rotation.w();
        T angle=orientation.norm();TV axis=(fabs(angle)>1e-8?orientation.normalized():TV::UnitX());
        T s=sin(angle/2);T s_over_angle=sinc(angle/2)/2,c=cos(angle/2);
    
        M_VxV axis_projection=axis*axis.transpose();
        M_VxV axis_orthogonal_projection=M_VxV::Identity()-axis_projection;
        M_VxV dudw=(c/2)*axis_projection+s_over_angle*axis_orthogonal_projection;
        TV dadw=-s/2*axis;
        return v*dadw.transpose()+cross_sign*b*dudw-Cross_Product_Matrix(v)*dudw;
    }

    static M_VxV Relative_Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,const ROTATION<TV>& target,TV& rotation_error_vector)
    {
        ROTATION<TV> composed_rotation=rotation*relative_rotation;
        rotation_error_vector=composed_rotation.vec()*composed_rotation.Sign()-target.Sign()*target.vec();
        return Compute_Orientation_Constraint_Matrix(rotation,relative_rotation,composed_rotation.Sign());
    }

    static M_VxV Orientation_Constraint_Matrix(const ROTATION<TV>& rotation,const ROTATION<TV>& relative_rotation,TV& rotation_error_vector)
    {
        ROTATION<TV> composed_rotation=rotation*relative_rotation;
        rotation_error_vector=composed_rotation.vec()*composed_rotation.Sign();
        return Compute_Orientation_Constraint_Matrix(rotation,relative_rotation,composed_rotation.Sign());
    }
};
}
#endif
