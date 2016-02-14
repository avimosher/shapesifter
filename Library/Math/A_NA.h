#ifndef __A_NA__
#define __A_NA__

#include <Math/ONE_NA.h>

namespace Mechanics{

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
}
#endif
