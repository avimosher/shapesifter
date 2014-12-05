//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class EVOLUTION
//#####################################################################
#ifndef __EVOLUTION__
#define __EVOLUTION__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class EVOLUTION_STEP;
template<class TV> class EVOLUTION_TYPE;
template<class TV> class FORCE;

template<class TV>
class EVOLUTION
{
    typedef typename TV::Scalar T;

    std::vector<EVOLUTION_TYPE<TV>*> evolution;
    std::vector<EVOLUTION_STEP<TV>*> steps;
public:
    EVOLUTION();
    ~EVOLUTION();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Evolution_Step(EVOLUTION_STEP<TV>& step,DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
    void Advance_One_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
};
}
#endif
