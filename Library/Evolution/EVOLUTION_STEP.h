#ifndef __EVOLUTION_STEP__
#define __EVOLUTION_STEP__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV>
class EVOLUTION_STEP
{
    typedef typename TV::Scalar T;

public:
    std::string name;
    T up_to_date_time;
    std::vector<EVOLUTION_STEP<TV>*> prerequisites;

    EVOLUTION_STEP():up_to_date_time(0){}
    ~EVOLUTION_STEP(){}

    bool Up_To_Date(SIMULATION<TV>& simulation,const T dt,const T time){
        return up_to_date_time>=time;
    }

    void Full_Step(SIMULATION<TV>& simulation,const T dt,const T time);
    virtual void Step(SIMULATION<TV>& simulation,const T dt,const T time)=0;
};
}
#endif
