//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class EQUATION
//#####################################################################
#ifndef __EQUATION__
#define __EQUATION__

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class EQUATION
{
    typedef typename TV::Scalar T;

    std::vector<EQUATION_TYPE*> evolution;
public:
    EQUATION();
    ~EQUATION();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Position_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
    void Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
};
}
#endif
