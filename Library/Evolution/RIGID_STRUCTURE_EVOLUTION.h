//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE_EVOLUTION
///////////////////////////////////////////////////////////////////////
#ifndef __RIGID_STRUCTURE_EVOLUTION__
#define __RIGID_STRUCTURE_EVOLUTION__

#include <Evolution/EVOLUTION_TYPE.h>

namespace Mechanics{
template<class TV> class RIGID_STRUCTURE;

template<class TV>
class RIGID_STRUCTURE_EVOLUTION:public EVOLUTION_TYPE<TV>
{
    typedef typename TV::Scalar T;

    std::vector<RIGID_STRUCTURE<TV>*> structures;
public:
    RIGID_STRUCTURE_EVOLUTION();
    ~RIGID_STRUCTURE_EVOLUTION();

    virtual T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T time,const T target_time);
    virtual void Position_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
    virtual void Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
};
}
#endif
