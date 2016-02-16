#ifndef __RXO__
#define __RXO__

#include <Math/Q_V.h>
#include <Math/Q_W.h>

namespace Mechanics{
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
        T norm_spin=std::max((T)epsilon(),spin[R].norm());
        TV_T dw_dspin=Q_W<TV>::First_Derivative(spin[R],norm_spin);
        M_VxV dq_dspin=Q_V<TV>::First_Derivative(spin[R],norm_spin);
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
}
#endif
