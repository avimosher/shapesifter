#ifndef __RCF_NF__
#define __RCF_NF__

#include <Math/F_NF.h>
#include <Math/RXO.h>

namespace Mechanics{
// r x f/|f|
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
}
#endif
