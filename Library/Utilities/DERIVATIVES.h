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

template<class TV,class FTYPE,class Derived>
struct Function
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;

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

        
    static TV Contract(const T_TENSOR& t,const TV& v1,const TV& v2){
        TV result;result.setZero();
        std::array<int,3> index{};
        for(index[0]=0;index[0]<3;index[0]++){
            for(index[1]=0;index[1]<3;index[1]++){
                for(index[2]=0;index[2]<3;index[2]++){
                    result(index[2])+=t(index[0],index[1],index[2])*v1(index[0])*v2(index[1]);}}}
        return result;
    }

    static T Contract(const M_VxV& m,const TV& v1,const TV& v2){
        return v1.transpose()*m*v2;
    }

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
        FTYPE initial=Derived::Evaluate(x,spin,offset);
        std::array<std::array<TV,2>,2> epsdx;
        std::array<TV,2> epsx;
        std::array<T_SPIN,2> epsspin;
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){epsdx[i][j]=eps*dx[i][j];}
            epsx[i]=x[i]+epsdx[i][0];
            epsspin[i]=spin[i]+epsdx[i][1];}
        FTYPE predicted=Apply_Derivatives(f,spin,offset,epsdx);
        FTYPE actual=Derived::Evaluate(epsx,epsspin,offset)-initial;
        return Norm(actual-predicted);
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

    template<int V1,int VTYPE,std::enable_if_t<VTYPE==LINEAR>* = nullptr>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return M_VxV::Identity()*VSIGN<V1>::SIGN;
    }

    template<int V1,int VTYPE,std::enable_if_t<VTYPE==ANGULAR>* = nullptr>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin[V1],offset[V1]).transpose()*VSIGN<V1>::SIGN;
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

    static T Evaluate(const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::Evaluate(x,spin,offset).norm();}

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

    static T Evaluate(const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return 1/F<TV>::Evaluate(x,spin,offset).norm();}

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

    static TV Evaluate(const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::Evaluate(x,spin,offset).normalized();}

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
        return BASE::Outer_Product(df_dv1,dnfinv_dv2,{0,2,1})+
            BASE::Outer_Product(df_dv2,dnfinv_dv1,{1,2,0})+
            BASE::Outer_Product(d2nfinv_dv2,f,{0,1,2})+
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

// rotation times offset
template<class TV,int R>
struct RXO:public Function<TV,TV,RXO<TV,R>>
{
    typedef Function<TV,TV,RXO<TV,R>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;

    static TV Evaluate(const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return ROTATION<TV>::From_Rotation_Vector(spin[R])*offset[R];}

    template<int V1,int VTYPE,std::enable_if_t<VTYPE==ANGULAR && V1==R>* = nullptr>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin[R],offset[R]).transpose();}

    template<int V1,int VTYPE,std::enable_if_t<!(VTYPE==ANGULAR && V1==R)>* = nullptr>
    static M_VxV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return M_VxV::Zero();}

    template<int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<VTYPE1==ANGULAR && VTYPE2==ANGULAR && V1==V2 && V1==R>* = nullptr>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return RIGID_STRUCTURE_INDEX_MAP<TV>::d2so_dSpin2(spin[R],offset[R]);
    }

    template<int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<!(VTYPE1==ANGULAR && VTYPE2==ANGULAR && V1==V2 && V1==R)>* = nullptr>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T_TENSOR t;t.setZero();return t;
    }
};

// R specifies whether the offset is V1 or V2
template<class TV,int R>
struct RCF_NF:public Function<TV,TV,RCF_NF<TV,R>>
{
    typedef Function<TV,TV,RCF_NF<TV,R>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;

    static TV Evaluate(const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return RXO<TV,R>::Evaluate(x,spin,offset).cross(F<TV>::Evaluate(x,spin,offset).normalized());}

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
        return BASE::Cross_Product(-Cross_Product_Matrix(f_nf),d2r_da2)+
            BASE::Cross_Product(dr_da,df_db)+
            BASE::Cross_Product(-df_da,dr_db)+
            BASE::Cross_Product(Cross_Product_Matrix(r),d2f_nf_da2);

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
