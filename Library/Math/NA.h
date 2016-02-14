#ifndef __NA__
#define __NA__

#include <Math/ONE_NA.h>

namespace Mechanics{

// |a|
template<class TV>
struct NA
{
    typedef typename TV::Scalar T;
    enum {d=TV::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;

    static M_VxV Second_Derivative(const TV& a,const T na){
        return M_VxV::Identity()/na+a*ONE_NA<TV>::First_Derivative(a,na).transpose();
    }
};
}
#endif
