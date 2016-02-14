#ifndef __NFINV__
#define __NFINV__

#include <Math/NF.h>

namespace Mechanics{
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
}
#endif
