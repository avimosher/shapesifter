#ifndef __F_NF__
#define __F_NF__

#include <Math/NFINV.h>

namespace Mechanics{
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
}
#endif
