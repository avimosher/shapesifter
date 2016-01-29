#ifndef __DERIVATIVES__
#define __DERIVATIVES__

#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <type_traits>

namespace std{
template<bool B,class T=void>
using enable_if_t=typename enable_if<B,T>::type;
}

namespace Mechanics{

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
    enum{LINEAR,ANGULAR};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;
        
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
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    enum{LINEAR,ANGULAR};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;

    F(){}

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
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    enum{LINEAR,ANGULAR};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;

    NF(){}

    static T Evaluate(const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::Evaluate(x,spin,offset).norm();
    }

    template<int V1,int VTYPE>
    static TV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::template First_Derivative<V1,VTYPE>(f,spin,offset)*f.normalized();
    }

    template<int V1,int VTYPE1,int V2,int VTYPE2>
        static M_VxV Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=f.norm();
        TV f_nf=f/nf;
        M_VxV df_dv1=F<TV>::template First_Derivative<V1,VTYPE1>(f,spin,offset).transpose();
        M_VxV df_dv2=F<TV>::template First_Derivative<V2,VTYPE2>(f,spin,offset).transpose();
        T_TENSOR d2f_dv2=F<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset);
        return 1/nf*df_dv1.transpose()*(M_VxV::Identity()-f_nf*f_nf.transpose())*df_dv2+
            RIGID_STRUCTURE_INDEX_MAP<TV>::Contract(d2f_dv2,f_nf,{1,2,0});
    }
};

}
#endif
