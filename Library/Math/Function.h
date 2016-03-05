#ifndef __Function__
#define __Function__

#include <Utilities/EIGEN_HELPERS.h>
#include <type_traits>

namespace std{
template<bool B,class T=void>
using enable_if_t=typename enable_if<B,T>::type;
}

namespace Mechanics{

enum VARTYPE{LINEAR,ANGULAR};
enum BODY{FIRST,SECOND};

#define WRAP_FUNCTION(function) function<
#define EXPAND_FUNCTION(function,...) function<__VA_ARGS__,
#define EXPAND_BODIES_TYPES(fixed,...) \
    fixed FIRST,LINEAR>(__VA_ARGS__);       \
    fixed FIRST,ANGULAR>(__VA_ARGS__);      \
    fixed SECOND,LINEAR>(__VA_ARGS__);      \
    fixed SECOND,ANGULAR>(__VA_ARGS__);

template<int V> struct VSIGN{};
template<> struct VSIGN<0>{enum{SIGN=-1};};
template<> struct VSIGN<1>{enum{SIGN=1};};

template<class TV> struct F;

template<class TV>
TV clamped_normalize(const TV& f){
    return f/std::max((typename TV::Scalar)epsilon(),f.norm());
}

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

    static T Norm(const TV& v){return v.norm();}
    static T Norm(const T v){return fabs(v);}

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
}
#endif
