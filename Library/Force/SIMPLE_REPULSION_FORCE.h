//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class SIMPLE_REPULSION_FORCE
//#####################################################################
#ifndef __SIMPLE_REPULSION_FORCE__
#define __SIMPLE_REPULSION_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class SIMPLE_REPULSION_FORCE
{
    typedef typename TV::Scalar T;

    std::vector<FORCE_TYPE*> evolution;
public:
    SIMPLE_REPULSION_FORCE();
    ~SIMPLE_REPULSION_FORCE();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Linearize(DATA<TV>& data,const T dt,const T target_time);
};
}
#endif
