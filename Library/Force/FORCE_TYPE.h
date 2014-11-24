//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class FORCE_TYPE
//#####################################################################
#ifndef __FORCE_TYPE__
#define __FORCE_TYPE__

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class FORCE_TYPE
{
    typedef typename TV::Scalar T;

public:
    FORCE_TYPE();
    ~FORCE_TYPE();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Position_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
    void Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
};
}
#endif
