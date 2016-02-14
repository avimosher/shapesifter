#ifndef __Q_W__
#define __Q_W__

#include <Math/NA.h>

namespace Mechanics{
// scalar part of quaternion from spin
template<class TV>
struct Q_W
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;

    static T Evaluate(const T_SPIN& spin,const T norm_spin){return cos(norm_spin/2);}
    static TV First_Derivative(const T_SPIN& spin,const T norm_spin){return -sinc(norm_spin/2)/4*spin.transpose();}

    static M_VxV Second_Derivative(const T_SPIN& spin,const T norm_spin){
        TV dna_da=spin/norm_spin;
        M_VxV d2na_da2=NA<TV>::Second_Derivative(spin,norm_spin);
        return -(T).25*cos(norm_spin/2)*dna_da*dna_da.transpose()-(T).5*sin(norm_spin/2)*d2na_da2;
    }
};
}
#endif
