//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class EVOLUTION_TYPE
///////////////////////////////////////////////////////////////////////
#ifndef __EVOLUTION_TYPE__
#define __EVOLUTION_TYPE__

#include <Data/EVOLUTION_TYPE.h>

namespace Mechanics{

template<class TV>
class EVOLUTION_TYPE
{
    typedef typename TV::Scalar T;

public:
    EVOLUTION_TYPE();
    ~EVOLUTION_TYPE();

    virtual void Position_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)=0;
    virtual void Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)=0;
    virtual void Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)=0;
};
}
#endif
