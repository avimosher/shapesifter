//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class EVOLUTION
//#####################################################################
#ifndef __EVOLUTION__
#define __EVOLUTION__

#include <Utilities/TYPE_UTILITIES.h>
#include <memory>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class EVOLUTION_STEP;
template<class TV> class EVOLUTION_TYPE;
template<class TV> class FORCE;

template<class TV>
class EVOLUTION : public std::vector<std::shared_ptr<EVOLUTION_STEP<TV>>>
{
    typedef typename TV::Scalar T;

public:
    EVOLUTION();
    ~EVOLUTION();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T time,const T target_time,bool& done);
    void Advance_One_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
};
}
#endif
