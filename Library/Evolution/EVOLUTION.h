//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class EVOLUTION
//#####################################################################
#ifndef __EVOLUTION__
#define __EVOLUTION__

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class EVOLUTION
{
    typedef typename TV::Scalar T;

    std::vector<EVOLUTION_TYPE*> evolution;
public:
    EVOLUTION();
    ~EVOLUTION();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Position_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
    void Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
};
}
#endif
