//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class NONINERTIAL_SYSTEM
//#####################################################################
#ifndef __NONINERTIAL_SYSTEM__
#define __NONINERTIAL_SYSTEM__

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class NONINERTIAL_SYSTEM
{
    typedef typename TV::Scalar T;

public:
    NONINERTIAL_SYSTEM();
    ~NONINERTIAL_SYSTEM();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Position_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
    void Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
};
}
#endif
