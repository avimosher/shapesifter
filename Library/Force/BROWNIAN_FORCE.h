//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class BROWNIAN_FORCE
//#####################################################################
#ifndef __BROWNIAN_FORCE__
#define __BROWNIAN_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class BROWNIAN_FORCE : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;

public:
    BROWNIAN_FORCE();
    ~BROWNIAN_FORCE();

    virtual T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    virtual void Add_Stochastic_Terms(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,VELOCITY<TV>& right_hand_side);
};
}
#endif
