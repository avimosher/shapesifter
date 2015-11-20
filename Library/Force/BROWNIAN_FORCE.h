#ifndef __BROWNIAN_FORCE__
#define __BROWNIAN_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class BROWNIAN_FORCE:public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    Matrix<T,Dynamic,1> stored_right_hand_side;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    T temperature;

    BROWNIAN_FORCE(){}
    ~BROWNIAN_FORCE(){}

    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    DEFINE_TYPE_NAME("BROWNIAN_FORCE")
};
}
#endif
