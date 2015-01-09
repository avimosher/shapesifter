//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class EVOLUTION_STEP
///////////////////////////////////////////////////////////////////////
#ifndef __EVOLUTION_STEP__
#define __EVOLUTION_STEP__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE;

template<class TV>
class EVOLUTION_STEP
{
    typedef typename TV::Scalar T;

public:
    T up_to_date_time;
    std::vector<EVOLUTION_STEP<TV>*> prerequisites;

    EVOLUTION_STEP();
    ~EVOLUTION_STEP();

    bool Up_To_Date(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time){
        return up_to_date_time>=time;
    }

    void Full_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
    virtual void Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)=0;
};
}
#endif
