#ifndef __F__
#define __F__

#include <Math/RXO.h>

namespace Mechanics{
// x_2-x_1
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
        //return M_VxV::Zero();
        return RXO<TV,V1>::template First_Derivative<V1,VTYPE>(f,spin,offset)*VSIGN<V1>::SIGN;
    }

    template<int V1,int VTYPE1,int V2,int VTYPE2>
    static T_TENSOR Second_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return RXO<TV,V1>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset)*(T)VSIGN<V1>::SIGN;
    }
};
}
#endif
