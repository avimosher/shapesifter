#ifndef __EVOLUTION_STEP__
#define __EVOLUTION_STEP__

#include <Utilities/TYPE_UTILITIES.h>
#include <json/json.h>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV>
class EVOLUTION_STEP
{
    typedef typename TV::Scalar T;

public:
    T up_to_date_time;
    std::vector<std::string> prerequisites;

    EVOLUTION_STEP():up_to_date_time(0){}
    virtual ~EVOLUTION_STEP(){}

    bool Up_To_Date(SIMULATION<TV>& simulation,const T dt,const T time){
        return up_to_date_time>=time;
    }

    virtual std::string Name() const
    {return "EVOLUTION_STEP";}

    void Full_Step(SIMULATION<TV>& simulation,const T dt,const T time);
    virtual void Step(SIMULATION<TV>& simulation,const T dt,const T time)=0;
    virtual void Finalize(Json::Value& node){}
    virtual bool Success(){return true;}
};
}
#endif
