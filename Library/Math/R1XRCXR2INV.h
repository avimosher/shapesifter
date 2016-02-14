#ifndef __R1XRCR2INV__
#define __R1XRCR2INV__

#include <Math/Q_V.h>
#include <Math/Q_W.h>

namespace Mechanics{
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
        T ns=std::max((T)epsilon(),spin[0].norm());
        TV dw_ds=Q_W<TV>::First_Derivative(spin[0],ns);
        M_VxV dq_ds=Q_V<TV>::First_Derivative(spin[0],ns);
        return dw_ds*vec.transpose()+right.w()*dq_ds.transpose()-(Cross_Product_Matrix(vec)*dq_ds).transpose();}

    template<int V,std::enable_if_t<V==1>* = nullptr>
    static M_VxV First_Derivative(const ROTATION<TV>& RC,const std::array<T_SPIN,2>& spin){
        ROTATION<TV> left=ROTATION<TV>::From_Rotation_Vector(spin[0])*RC;
        ROTATION<TV> right=ROTATION<TV>::From_Rotation_Vector(spin[1]);
        TV vec=left.vec();
        T ns=std::max((T)epsilon(),spin[1].norm());
        TV dw_ds=Q_W<TV>::First_Derivative(spin[1],ns);
        M_VxV dq_ds=Q_V<TV>::First_Derivative(spin[1],ns);
        return dw_ds*vec.transpose()-left.w()*dq_ds.transpose()-(Cross_Product_Matrix(vec)*dq_ds).transpose();}

    static void Constraint_Derivatives(const int constraint_index,const std::array<int,2>& index,const ROTATION<TV>& RC,const std::array<T_SPIN,2>& spin,std::vector<Triplet<T>>& constraint_terms)
    {
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(constraint_index,index[0],1,1,First_Derivative<0>(RC,spin).transpose(),constraint_terms);
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(constraint_index,index[1],1,1,First_Derivative<1>(RC,spin).transpose(),constraint_terms);
    }
};
}
#endif
