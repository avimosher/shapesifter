#ifndef __Q_V__
#define __Q_V__

#include <Math/A_NA.h>

namespace Mechanics{
// vector part of quaternion from spin
template<class TV>
struct Q_V
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<d,d,d>> T_TENSOR;

    static TV Evaluate(const T_SPIN& spin,const T norm_spin){return sinc(norm_spin/2)*spin/2;}

    static M_VxV First_Derivative(const T_SPIN& spin,const T norm_spin){
        T_SPIN spin_normspin=spin/norm_spin;
        return cos(norm_spin/2)/2*spin_normspin*spin_normspin.transpose()+sinc(norm_spin/2)/2*(Matrix<T,t,t>::Identity()-spin_normspin*spin_normspin.transpose());
    }

    static T_TENSOR Second_Derivative(const T_SPIN& spin,const T norm_spin){
        T_SPIN s_ns=spin/norm_spin;
        T_SPIN dns_ds=s_ns;
        M_VxV ds_ns_ds=A_NA<TV>::First_Derivative(spin,norm_spin);
        M_VxV d2ns_ds2=ds_ns_ds;
        T_TENSOR d2s_ns_ds2=A_NA<TV>::Second_Derivative(spin,norm_spin);
        
        T s=sinc(norm_spin/2);
        T c=cos(norm_spin/2);
        return Outer_Product(dns_ds*dns_ds.transpose(),spin,{0,1,2})*(-(T).125*s)+
            (Outer_Product(ds_ns_ds,dns_ds,{2,0,1})+Outer_Product(ds_ns_ds,dns_ds,{2,1,0})+Outer_Product(d2ns_ds2,s_ns,{0,1,2}))*(T).5*c+
            d2s_ns_ds2*(s*norm_spin/2);
    }
};
}
#endif
