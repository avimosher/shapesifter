#include <Evolution/EVOLUTION_STEP.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void EVOLUTION_STEP<TV>::
Full_Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    for(auto prerequisite_step : prerequisites) {
        if(!prerequisite_step->Up_To_Date(simulation,dt,time)) {
            prerequisite_step->Full_Step(simulation,dt,time);}
    }
    if(!Up_To_Date(simulation,dt,time)){
        Step(simulation,dt,time);
        up_to_date_time=time;}
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EVOLUTION_STEP)
