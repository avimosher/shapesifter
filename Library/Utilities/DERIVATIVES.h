#ifndef __DERIVATIVES__
#define __DERIVATIVES__

#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <type_traits>

namespace std{
template<bool B,class T=void>
using enable_if_t=typename enable_if<B,T>::type;
}

namespace Mechanics{

enum VARTYPE{LINEAR,ANGULAR};

template<int V> struct VSIGN{};
template<> struct VSIGN<0>{enum{SIGN=-1};};
template<> struct VSIGN<1>{enum{SIGN=1};};

template<class TV> struct F;

inline double epsilon(){return 1e-8;}

template<class TV,class FTYPE,class Derived>
struct Function
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,1,d> TV_T;
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;

    template<int V1,int VTYPE1,int V2,int VTYPE2>
    static FTYPE Apply_Second(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<std::array<TV,2>,2>& dx)
    {
        return Contract(Derived::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),dx[V1][VTYPE1],dx[V2][VTYPE2]);
    }

    template<int V1,int VTYPE1>
    static FTYPE Apply_First(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<std::array<TV,2>,2>& dx)
    {
        return Derived::template First_Derivative<V1,VTYPE1>(f,spin,offset).transpose()*dx[V1][VTYPE1];
    }

    template<int V1,int VTYPE1>
    static FTYPE Apply_Part_Two(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<std::array<TV,2>,2>& dx)
    {
        return Apply_Second<V1,VTYPE1,0,LINEAR>(f,spin,offset,dx)+Apply_Second<V1,VTYPE1,0,ANGULAR>(f,spin,offset,dx)+
            Apply_Second<V1,VTYPE1,1,LINEAR>(f,spin,offset,dx)+Apply_Second<V1,VTYPE1,1,ANGULAR>(f,spin,offset,dx);
    }

    template<int V1,int VTYPE1>
    static FTYPE Apply_Part_One(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<std::array<TV,2>,2>& dx)
    {
        return Apply_First<V1,VTYPE1>(f,spin,offset,dx)+(T).5*Apply_Part_Two<V1,VTYPE1>(f,spin,offset,dx);
    }

    static FTYPE Apply_Derivatives(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<std::array<TV,2>,2>& dx)
    {
        return Apply_Part_One<0,LINEAR>(f,spin,offset,dx)+Apply_Part_One<0,ANGULAR>(f,spin,offset,dx)+
            Apply_Part_One<1,LINEAR>(f,spin,offset,dx)+Apply_Part_One<1,ANGULAR>(f,spin,offset,dx);
    }

    static T Norm(const TV& v){
        return v.norm();
    }

    static T Norm(const T v){
        return fabs(v);
    }

    static T Test_Error(const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<std::array<TV,2>,2>& dx,T eps){
        TV f=F<TV>::Evaluate(x,spin,offset);
        FTYPE initial=Derived::Evaluate(f,spin,offset);
        std::array<std::array<TV,2>,2> epsdx;
        std::array<TV,2> epsx;
        std::array<T_SPIN,2> epsspin;
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){epsdx[i][j]=eps*dx[i][j];}
            epsx[i]=x[i]+epsdx[i][0];
            epsspin[i]=spin[i]+epsdx[i][1];}
        FTYPE predicted=Apply_Derivatives(f,spin,offset,epsdx);
        TV f_final=F<TV>::Evaluate(epsx,epsspin,offset);
        FTYPE actual=Derived::Evaluate(f_final,epsspin,offset)-initial;
        return Norm(actual-predicted);
    }
};



// 1/|a|
template<class TV>
struct ONE_NA
{
    typedef typename TV::Scalar T;
    static T Evaluate(const TV& a,const T na){return 1/na;}
    static TV First_Derivative(const TV& a,const T na){return -a/cube(na);}
};


// |a|
template<class TV>
struct NA
{
    typedef typename TV::Scalar T;
    enum {d=TV::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;

    static M_VxV Second_Derivative(const TV& a,const T na){
        return M_VxV::Identity()/na+a*ONE_NA<TV>::First_Derivative(a,na).transpose();
    }
};

// a/|a|
template<class TV>
struct A_NA
{
    typedef typename TV::Scalar T;
    enum {d=TV::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<d,d,d>> T_TENSOR;

    static TV Evaluate(const TV& a,const T na){return a/na;}

    static M_VxV First_Derivative(const TV& a,const T na){
        if(na<epsilon()){return M_VxV::Zero();}
        T one_na=1/na;
        TV a_na=a*one_na;
        return one_na*(M_VxV::Identity()-a_na*a_na.transpose());
    }

    static T_TENSOR Second_Derivative(const TV& a,const T na){
        M_VxV da_na_da=First_Derivative(a,na);
        T one_na=1/na;
        TV a_na=a*one_na;
        TV dnainv_da=ONE_NA<TV>::First_Derivative(a,na);
        
        return Outer_Product(M_VxV::Identity()-a_na*a_na.transpose(),dnainv_da,{2,0,1})-
            Outer_Product(da_na_da,a_na,{2,0,1})*(1/na)-Outer_Product(da_na_da,a_na,{1,0,2})*(1/na);
    }
};

// scalar part of quaternion from spin
template<class TV>
struct Q_W
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;

    static T Evaluate(const T_SPIN& spin,const T norm_spin){return cos(norm_spin/2);}
    static TV First_Derivative(const T_SPIN& spin,const T norm_spin){return -sinc(norm_spin/2)/4*spin.transpose();}

    static M_VxV Second_Derivative(const T_SPIN& spin,const T norm_spin){
        TV dna_da=spin/norm_spin;
        M_VxV d2na_da2=NA<TV>::Second_Derivative(spin,norm_spin);
        return -(T).25*cos(norm_spin/2)*dna_da*dna_da.transpose()-(T).5*sin(norm_spin/2)*d2na_da2;
    }
};

// vector part of quaternion from spin
template<class TV>
struct Q_V
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<d,d,d>> T_TENSOR;

    static TV Evaluate(const T_SPIN& spin,const T norm_spin){return sinc(norm_spin/2)*spin/2;}

    static M_VxV First_Derivative(const T_SPIN& spin,const T norm_spin){
        T_SPIN spin_normspin=spin/norm_spin;
        return cos(norm_spin/2)/2*spin_normspin*spin_normspin.transpose()+sinc(norm_spin/2)/2*(Matrix<T,t,t>::Identity()-spin_normspin*spin_normspin.transpose());
    }

    static T_TENSOR Second_Derivative(const T_SPIN& spin,const T norm_spin){
        T_SPIN s_ns=spin/norm_spin;
        T_SPIN dns_ds=s_ns;
        M_VxV ds_ns_ds=A_NA<TV>::First_Derivative(spin,norm_spin);
        M_VxV d2ns_ds2=ds_ns_ds;
        T_TENSOR d2s_ns_ds2=A_NA<TV>::Second_Derivative(spin,norm_spin);
        
        T s=sinc(norm_spin/2);
        T c=cos(norm_spin/2);
        return Outer_Product(dns_ds*dns_ds.transpose(),spin,{0,1,2})*(-(T).125*s)+
            (Outer_Product(ds_ns_ds,dns_ds,{2,0,1})+Outer_Product(ds_ns_ds,dns_ds,{2,1,0})+Outer_Product(d2ns_ds2,s_ns,{0,1,2}))*(T).5*c+
            d2s_ns_ds2*(s*norm_spin/2);
    }
};

// rotation times offset
template<class TV,int R>
struct RXO:public Function<TV,TV,RXO<TV,R>>
{
    typedef Function<TV,TV,RXO<TV,R>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;using typename BASE::TV_T;

    static TV Evaluate(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return ROTATION<TV>::From_Rotation_Vector(spin[R])*offset[R];}

    template<int V1,int VTYPE,std::enable_if_t<VTYPE==ANGULAR && V1==R>* = nullptr>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T norm_spin=spin[R].norm();
        TV_T dw_dspin=Q_W<TV>::First_Derivative(spin[R],norm_spin);
        TV spin_normspin=(norm_spin>epsilon())?(TV)(spin[R]/norm_spin):TV::UnitX();
        M_VxV dq_dspin=Q_V<TV>::First_Derivative(spin_normspin,norm_spin);
        T w=cos(norm_spin/2);
        TV q=sinc(norm_spin/2)*spin[R]/2;
        return (2*q.cross(offset[R])*dw_dspin-2*(Cross_Product_Matrix(offset[R])*w+Cross_Product_Matrix(q.cross(offset[R]))+Cross_Product_Matrix(q)*Cross_Product_Matrix(offset[R]))*dq_dspin).transpose();
    }

    template<int V1,int VTYPE,std::enable_if_t<!(VTYPE==ANGULAR && V1==R)>* = nullptr>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return M_VxV::Zero();}

    template<int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<VTYPE1==ANGULAR && VTYPE2==ANGULAR && V1==V2 && V1==R>* = nullptr>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        ROTATION<TV> rotation(ROTATION<TV>::From_Rotation_Vector(spin[R]));
        TV q=rotation.vec();
        T w=rotation.w();
        M_VxV ostar=Cross_Product_Matrix(offset[R]);
        T ns=std::max((T)epsilon(),spin[R].norm());
        M_VxV d2w_ds2=Q_W<TV>::Second_Derivative(spin[R],ns);
        T_TENSOR d2q_ds2=Q_V<TV>::Second_Derivative(spin[R],ns);
        TV s_ns=spin[R]/ns;
        M_VxV dq_ds=Q_V<TV>::First_Derivative(spin[R],ns);
        TV dw_ds=Q_W<TV>::First_Derivative(spin[R],ns);
        return Outer_Product(-2*ostar*dq_ds,dw_ds,{2,1,0})+
            Outer_Product(d2w_ds2,(2*q.cross(offset[R])).eval(),{0,1,2})+
            Outer_Product(ostar*dq_ds,(dw_ds*(-2)).eval(),{2,0,1})+Cross_Product(-2*dq_ds,(ostar*dq_ds).eval())+Cross_Product(2*ostar*dq_ds,dq_ds)+
            Cross_Product(-2*(w*ostar+Cross_Product_Matrix(q.cross(offset[R]))+Cross_Product_Matrix(q)*ostar),d2q_ds2);
    }

    template<int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<!(VTYPE1==ANGULAR && VTYPE2==ANGULAR && V1==V2 && V1==R)>* = nullptr>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T_TENSOR t;t.setZero();return t;
    }
};


template<class TV>
struct F:public Function<TV,TV,F<TV>>
{
    typedef Function<TV,TV,F<TV>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;

    static TV Evaluate(const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return x[1]+ROTATION<TV>::From_Rotation_Vector(spin[1])*offset[1]-(x[0]+ROTATION<TV>::From_Rotation_Vector(spin[0])*offset[0]);
    }

    static TV Evaluate(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return f;
    }

    template<int V1,int VTYPE,std::enable_if_t<VTYPE==LINEAR>* = nullptr>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return M_VxV::Identity()*VSIGN<V1>::SIGN;
    }

    template<int V1,int VTYPE,std::enable_if_t<VTYPE==ANGULAR>* = nullptr>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return RXO<TV,V1>::template First_Derivative<V1,VTYPE>(f,spin,offset)*VSIGN<V1>::SIGN;
    }

    template<int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<VTYPE1==ANGULAR && VTYPE2==ANGULAR && V1==V2>* = nullptr>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return RIGID_STRUCTURE_INDEX_MAP<TV>::d2so_dSpin2(spin[V1],offset[V1])*(T)VSIGN<V1>::SIGN;
    }

    template<int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<!(VTYPE1==ANGULAR && VTYPE2==ANGULAR && V1==V2)>* = nullptr>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T_TENSOR t;t.setZero();return t;
    }
};

template<class TV>
struct NF:public Function<TV,typename TV::Scalar,NF<TV>>
{
    typedef Function<TV,typename TV::Scalar,NF<TV>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;

    static T Evaluate(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::Evaluate(f,spin,offset).norm();}

    template<int V1,int VTYPE>
    static TV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::template First_Derivative<V1,VTYPE>(f,spin,offset)*f.normalized();}

    template<int V1,int VTYPE1,int V2,int VTYPE2>
        static M_VxV Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=std::max(f.norm(),(T)1e-8);
        TV f_nf=f/nf;
        M_VxV df_dv1=F<TV>::template First_Derivative<V1,VTYPE1>(f,spin,offset).transpose();
        M_VxV df_dv2=F<TV>::template First_Derivative<V2,VTYPE2>(f,spin,offset).transpose();
        T_TENSOR d2f_dv2=F<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset);
        return 1/nf*df_dv1.transpose()*(M_VxV::Identity()-f_nf*f_nf.transpose())*df_dv2+
            RIGID_STRUCTURE_INDEX_MAP<TV>::Contract(d2f_dv2,f_nf,{1,2,0});
    }
};

template<class TV>
struct NFINV:public Function<TV,typename TV::Scalar,NFINV<TV>>
{
    typedef Function<TV,typename TV::Scalar,NFINV<TV>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;

    static T Evaluate(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return 1/F<TV>::Evaluate(f,spin,offset).norm();}

    template<int V1,int VTYPE>
    static TV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return -F<TV>::template First_Derivative<V1,VTYPE>(f,spin,offset)*f/cube(std::max((T)1e-8,f.norm()));}

    template<int V1,int VTYPE1,int V2,int VTYPE2>
        static M_VxV Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=std::max((T)1e-8,f.norm());
        TV dnf_dv1=NF<TV>::template First_Derivative<V1,VTYPE1>(f,spin,offset);
        TV dnf_dv2=NF<TV>::template First_Derivative<V2,VTYPE2>(f,spin,offset);
        M_VxV d2nf_dv2=NF<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset);
        return 2/cube(nf)*dnf_dv1*dnf_dv2.transpose()-1/sqr(nf)*d2nf_dv2;
    }
};

template<class TV>
struct F_NF:public Function<TV,TV,F_NF<TV>>
{
    typedef Function<TV,TV,F_NF<TV>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;

    static TV Evaluate(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::Evaluate(f,spin,offset).normalized();}

    template<int V1,int VTYPE>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::template First_Derivative<V1,VTYPE>(f,spin,offset)/std::max((T)1e-8,f.norm())+NFINV<TV>::template First_Derivative<V1,VTYPE>(f,spin,offset)*f.transpose();}

    template<int V1,int VTYPE1,int V2,int VTYPE2>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=std::max((T)1e-8,f.norm());
        M_VxV df_dv1=F<TV>::template First_Derivative<V1,VTYPE1>(f,spin,offset);
        M_VxV df_dv2=F<TV>::template First_Derivative<V2,VTYPE2>(f,spin,offset);
        TV dnfinv_dv1=NFINV<TV>::template First_Derivative<V1,VTYPE1>(f,spin,offset);
        TV dnfinv_dv2=NFINV<TV>::template First_Derivative<V2,VTYPE2>(f,spin,offset);
        M_VxV d2nfinv_dv2=NFINV<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset);
        T_TENSOR d2f_dv2=F<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset);
        return Outer_Product(df_dv1,dnfinv_dv2,{0,2,1})+
            Outer_Product(df_dv2,dnfinv_dv1,{1,2,0})+
            Outer_Product(d2nfinv_dv2,f,{0,1,2})+
            d2f_dv2*(1/nf);
    }
};

template<class TV>
struct R1XRCXR2INV
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;

    static TV Evaluate(const ROTATION<TV>& RC,const std::array<T_SPIN,2>& spin){
        return (ROTATION<TV>::From_Rotation_Vector(spin[0])*RC*ROTATION<TV>::From_Rotation_Vector(spin[1]).inverse()).vec();}

    template<int V,std::enable_if_t<V==0>* = nullptr>
    static M_VxV First_Derivative(const ROTATION<TV>& RC,const std::array<T_SPIN,2>& spin){
        ROTATION<TV> right=RC*ROTATION<TV>::From_Rotation_Vector(spin[1]).inverse();
        ROTATION<TV> left=ROTATION<TV>::From_Rotation_Vector(spin[0]);
        TV vec=right.vec();
        T ns=std::max((T)1e-8,spin[0].norm());
        TV dw_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::dw_dSpin(spin[0],ns);
        M_VxV dq_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::dq_dSpin(spin[0]/ns,ns);
        return dw_ds*vec.transpose()+right.w()*dq_ds.transpose()-(Cross_Product_Matrix(vec)*dq_ds).transpose();}

    template<int V,std::enable_if_t<V==1>* = nullptr>
    static M_VxV First_Derivative(const ROTATION<TV>& RC,const std::array<T_SPIN,2>& spin){
        ROTATION<TV> left=ROTATION<TV>::From_Rotation_Vector(spin[0])*RC;
        ROTATION<TV> right=ROTATION<TV>::From_Rotation_Vector(spin[1]);
        TV vec=left.vec();
        T ns=std::max((T)1e-8,spin[1].norm());
        TV dw_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::dw_dSpin(spin[1],ns);
        M_VxV dq_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::dq_dSpin(spin[1]/ns,ns);
        return dw_ds*vec.transpose()-left.w()*dq_ds.transpose()-(Cross_Product_Matrix(vec)*dq_ds).transpose();}

    static void Constraint_Derivatives(const int constraint_index,const std::array<int,2>& index,const ROTATION<TV>& RC,const std::array<T_SPIN,2>& spin,std::vector<Triplet<T>>& constraint_terms)
    {
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(constraint_index,index[0],1,1,First_Derivative<0>(RC,spin).transpose(),constraint_terms);
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(constraint_index,index[1],1,1,First_Derivative<1>(RC,spin).transpose(),constraint_terms);
    }
};

// R specifies whether the offset is V1 or V2
template<class TV,int R>
struct RCF_NF:public Function<TV,TV,RCF_NF<TV,R>>
{
    typedef Function<TV,TV,RCF_NF<TV,R>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;

    static TV Evaluate(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return RXO<TV,R>::Evaluate(f,spin,offset).cross(F<TV>::Evaluate(f,spin,offset).normalized());}

    // this is only valid when V1 and R match and VTYPE is ANGULAR.  Otherwise the dr_da part is zero.
    template<int V1,int VTYPE>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        M_VxV dr_da=RXO<TV,R>::template First_Derivative<V1,VTYPE>(f,spin,offset).transpose();
        M_VxV df_da=F_NF<TV>::template First_Derivative<V1,VTYPE>(f,spin,offset).transpose();
        return (-Cross_Product_Matrix(f.normalized())*dr_da+Cross_Product_Matrix(ROTATION<TV>::From_Rotation_Vector(spin[R])*offset[R])*df_da).transpose();}

    template<int V1,int VTYPE1,int V2,int VTYPE2>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        M_VxV dr_da=RXO<TV,R>::template First_Derivative<V1,VTYPE1>(f,spin,offset).transpose();
        M_VxV dr_db=RXO<TV,R>::template First_Derivative<V2,VTYPE2>(f,spin,offset).transpose();
        M_VxV df_da=F_NF<TV>::template First_Derivative<V1,VTYPE1>(f,spin,offset).transpose();
        M_VxV df_db=F_NF<TV>::template First_Derivative<V2,VTYPE2>(f,spin,offset).transpose();
        TV f_nf=f.normalized();
        TV r=ROTATION<TV>::From_Rotation_Vector(spin[R])*offset[R];
        T_TENSOR d2r_da2=RXO<TV,R>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset);
        T_TENSOR d2f_nf_da2=F_NF<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset);
        return Cross_Product(-Cross_Product_Matrix(f_nf),d2r_da2)+
            Cross_Product(dr_da,df_db)+
            Cross_Product(-df_da,dr_db)+
            Cross_Product(Cross_Product_Matrix(r),d2f_nf_da2);
    }
};


template<class TV>
struct Spring_Force
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;

    template<int E,int ETYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static TV Evaluate(const T k,const T target,const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        TV f=F<TV>::Evaluate(x,spin,offset);
        T nf=std::max((T)1e-8,f.norm());
        return -k*(nf-target)*f/nf*VSIGN<E>::SIGN;
    }

    template<int E,int ETYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static T_SPIN Evaluate(const T k,const T target,const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        TV f=F<TV>::Evaluate(x,spin,offset);
        T nf=std::max((T)1e-8,f.norm());
        return (ROTATION<TV>::From_Rotation_Vector(spin[E])*offset[E]).cross(-k*(nf-target)*f/nf*VSIGN<E>::SIGN);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static M_VxV First_Derivative(const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=std::max((T)1e-8,f.norm());
        return -k*VSIGN<E>::SIGN*(NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset)*F_NF<TV>::Evaluate(f,spin,offset).transpose()+(nf-target)*F_NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset));
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static M_VxV First_Derivative(const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=std::max((T)1e-8,f.norm());
        return -k*VSIGN<E>::SIGN*(RCF_NF<TV,E>::template First_Derivative<V,VTYPE>(f,spin,offset)*(nf-target)+NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset)*RCF_NF<TV,E>::Evaluate(f,spin,offset).transpose());
    }

    template<int E,int ETYPE,int V,int VTYPE>
    static void Matrix_Term(const std::array<int,2>& index,const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[E],index[V],ETYPE,VTYPE,First_Derivative<E,ETYPE,V,VTYPE>(k,target,f,spin,offset).transpose(),force_terms);
    }

    template<int E,int ETYPE>
    static void First_Derivatives(const std::array<int,2>& index,const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Matrix_Term<E,ETYPE,0,LINEAR>(index,k,target,f,spin,offset,force_terms);
        Matrix_Term<E,ETYPE,0,ANGULAR>(index,k,target,f,spin,offset,force_terms);
        Matrix_Term<E,ETYPE,1,LINEAR>(index,k,target,f,spin,offset,force_terms);
        Matrix_Term<E,ETYPE,1,ANGULAR>(index,k,target,f,spin,offset,force_terms);
    }

    static void Derivatives(const std::array<int,2>& index,const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        First_Derivatives<0,LINEAR>(index,k,target,f,spin,offset,force_terms);
        First_Derivatives<0,ANGULAR>(index,k,target,f,spin,offset,force_terms);
        First_Derivatives<1,LINEAR>(index,k,target,f,spin,offset,force_terms);
        First_Derivatives<1,ANGULAR>(index,k,target,f,spin,offset,force_terms);
    }

    template<int E,int ETYPE>
    static void Evaluate_Term(const std::array<int,2>& index,const T k,const T target,const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& error){
        error.template block<d,1>(index[E]*(t+d)+ETYPE*d,0)+=Evaluate<E,ETYPE>(k,target,x,spin,offset);
    }

    static void Evaluate(const std::array<int,2>& index,const T k,const T target,const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& error){
        Evaluate_Term<0,LINEAR>(index,k,target,x,spin,offset,error);
        Evaluate_Term<0,ANGULAR>(index,k,target,x,spin,offset,error);
        Evaluate_Term<1,LINEAR>(index,k,target,x,spin,offset,error);
        Evaluate_Term<1,ANGULAR>(index,k,target,x,spin,offset,error);
    }
};


template<class TV>
struct Relative_Position_Force
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;

    static M_VxV Contract(const T_TENSOR& t,const TV& v){
        M_VxV result;result.setZero();
        std::array<int,3> index{};
        for(index[0]=0;index[0]<3;index[0]++){
            for(index[1]=0;index[1]<3;index[1]++){
                for(index[2]=0;index[2]<3;index[2]++){
                    result(index[0],index[1])+=t(index[0],index[1],index[2])*v(index[2]);}}}
        return result;
    }

    template<int E,int ETYPE,int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void Apply(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms)
    {
        // do the F_NF version
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,Contract(F_NF<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),VSIGN<E>::SIGN*error[E][ETYPE]),hessian_terms);
    }

    template<int E,int ETYPE,int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void Apply(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms)
    {
        // do the RCF_NF version
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,Contract(RCF_NF<TV,E>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),VSIGN<E>::SIGN*error[E][ETYPE]),hessian_terms);
    }

    template<int E,int ETYPE,int V1,int VTYPE1>
    static void Apply_Part_One(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        Apply<E,ETYPE,V1,VTYPE1,0,LINEAR>(f,spin,offset,index,error,hessian_terms);
        Apply<E,ETYPE,V1,VTYPE1,0,ANGULAR>(f,spin,offset,index,error,hessian_terms);
        Apply<E,ETYPE,V1,VTYPE1,1,LINEAR>(f,spin,offset,index,error,hessian_terms);
        Apply<E,ETYPE,V1,VTYPE1,1,ANGULAR>(f,spin,offset,index,error,hessian_terms);
    }

    template<int E,int ETYPE>
    static void Derivatives(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        Apply_Part_One<E,ETYPE,0,LINEAR>(f,spin,offset,index,error,hessian_terms);
        Apply_Part_One<E,ETYPE,0,ANGULAR>(f,spin,offset,index,error,hessian_terms);
        Apply_Part_One<E,ETYPE,1,LINEAR>(f,spin,offset,index,error,hessian_terms);
        Apply_Part_One<E,ETYPE,1,ANGULAR>(f,spin,offset,index,error,hessian_terms);
    }

    static void Build(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        // do both equation sets
        Derivatives<0,LINEAR>(f,spin,offset,index,error,hessian_terms);
        Derivatives<0,ANGULAR>(f,spin,offset,index,error,hessian_terms);
        Derivatives<1,LINEAR>(f,spin,offset,index,error,hessian_terms);
        Derivatives<1,ANGULAR>(f,spin,offset,index,error,hessian_terms);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void Force_Velocity_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        TV term=(F_NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset))*error[E][ETYPE]*(T)VSIGN<E>::SIGN;
        Flatten_Matrix_Term<T,t+d,1,d,1>(index[V],constraint_index,VTYPE,0,term,force_constraint_terms);
        Flatten_Matrix_Term<T,1,t+d,1,d>(constraint_index,index[V],0,VTYPE,term.transpose(),constraint_force_terms);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void Force_Velocity_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        TV term=(RCF_NF<TV,E>::template First_Derivative<V,VTYPE>(f,spin,offset))*error[E][ETYPE]*(T)VSIGN<E>::SIGN;
        Flatten_Matrix_Term<T,t+d,1,d,1>(index[V],constraint_index,VTYPE,0,term,force_constraint_terms);
        Flatten_Matrix_Term<T,1,t+d,1,d>(constraint_index,index[V],0,VTYPE,term.transpose(),constraint_force_terms);
    }

    template<int E,int ETYPE>
    static void Force_Velocity_Derivatives(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        Force_Velocity_Derivative<E,ETYPE,0,LINEAR>(f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
        Force_Velocity_Derivative<E,ETYPE,0,ANGULAR>(f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
        Force_Velocity_Derivative<E,ETYPE,1,LINEAR>(f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
        Force_Velocity_Derivative<E,ETYPE,1,ANGULAR>(f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
    }

    static void Force_Velocity(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        Force_Velocity_Derivatives<0,LINEAR>(f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
        Force_Velocity_Derivatives<0,ANGULAR>(f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
        Force_Velocity_Derivatives<1,LINEAR>(f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
        Force_Velocity_Derivatives<1,ANGULAR>(f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void First_Derivative(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[E],index[V],ETYPE,VTYPE,F_NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset).transpose()*scale*(T)VSIGN<E>::SIGN,force_terms);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void First_Derivative(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[E],index[V],ETYPE,VTYPE,RCF_NF<TV,E>::template First_Derivative<V,VTYPE>(f,spin,offset).transpose()*scale*(T)VSIGN<E>::SIGN,force_terms);
    }

    template<int E,int ETYPE>
    static void First_Derivatives_Second_Part(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        First_Derivative<E,ETYPE,0,LINEAR>(index,scale,f,spin,offset,force_terms);
        First_Derivative<E,ETYPE,0,ANGULAR>(index,scale,f,spin,offset,force_terms);
        First_Derivative<E,ETYPE,1,LINEAR>(index,scale,f,spin,offset,force_terms);
        First_Derivative<E,ETYPE,1,ANGULAR>(index,scale,f,spin,offset,force_terms);
    }

    static void First_Derivatives(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        First_Derivatives_Second_Part<0,LINEAR>(index,scale,f,spin,offset,force_terms);
        First_Derivatives_Second_Part<0,ANGULAR>(index,scale,f,spin,offset,force_terms);
        First_Derivatives_Second_Part<1,LINEAR>(index,scale,f,spin,offset,force_terms);
        First_Derivatives_Second_Part<1,ANGULAR>(index,scale,f,spin,offset,force_terms);
    }

    template<int V1,int VTYPE1,int V2,int VTYPE2>
    static void Constraint_Second_Derivative(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,error*NF<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),hessian_terms);
    }

    template<int V1,int VTYPE1>
    static void Constraint_Second_Derivative_Part_Two(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        Constraint_Second_Derivative<V1,VTYPE1,0,LINEAR>(index,error,f,spin,offset,hessian_terms);
        Constraint_Second_Derivative<V1,VTYPE1,0,ANGULAR>(index,error,f,spin,offset,hessian_terms);
        Constraint_Second_Derivative<V1,VTYPE1,1,LINEAR>(index,error,f,spin,offset,hessian_terms);
        Constraint_Second_Derivative<V1,VTYPE1,1,ANGULAR>(index,error,f,spin,offset,hessian_terms);
    }
    
    static void Constraint_Second_Derivatives(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        Constraint_Second_Derivative_Part_Two<0,LINEAR>(index,error,f,spin,offset,hessian_terms);
        Constraint_Second_Derivative_Part_Two<0,ANGULAR>(index,error,f,spin,offset,hessian_terms);
        Constraint_Second_Derivative_Part_Two<1,LINEAR>(index,error,f,spin,offset,hessian_terms);
        Constraint_Second_Derivative_Part_Two<1,ANGULAR>(index,error,f,spin,offset,hessian_terms);
    }
};



}
#endif
