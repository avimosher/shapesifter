#ifndef __NF__
#define __NF__

#include <Math/F.h>

namespace Mechanics{
template<class TV>
struct NF:public Function<TV,typename TV::Scalar,NF<TV>>
{
    typedef Function<TV,typename TV::Scalar,NF<TV>> BASE;
    using typename BASE::T;using typename BASE::M_VxV;using typename BASE::T_TENSOR;using typename BASE::T_SPIN;

    static T Evaluate(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::Evaluate(f,spin,offset).norm();}

    template<int V1,int VTYPE>
    static TV First_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return F<TV>::template First_Derivative<V1,VTYPE>(f,spin,offset)*clamped_normalize(f);}

    template<int V1,int VTYPE1,int V2,int VTYPE2>
    static M_VxV Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=std::max(f.norm(),(T)1e-8);
        TV f_nf=f/nf;
        M_VxV df_dv1=F<TV>::template First_Derivative<V1,VTYPE1>(f,spin,offset).transpose();
        M_VxV df_dv2=F<TV>::template First_Derivative<V2,VTYPE2>(f,spin,offset).transpose();
        T_TENSOR d2f_dv2=F<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset);
        return df_dv1.transpose()*(M_VxV::Identity()-f_nf*f_nf.transpose())*df_dv2*(1/nf)+
            Contract(d2f_dv2,f_nf);
    }
};
}
#endif
