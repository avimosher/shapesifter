#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/EQUATION.h>
#include <Evolution/EVOLUTION.h>
#include <Evolution/EVOLUTION_STEP.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar EVOLUTION<TV>::
Compute_Dt(SIMULATION<TV>& simulation,const T time,const T target_time,bool& done)
{
    /*TIME_STEP_INFLUENCE<TV> time_step;
    done=false;
    simulation.force.Compute_Dt(time,simulation.data,time_step);
    T dt=max(simulation.dt,time_step.Compute_Dt(minimum_dt));
    if(time+dt>=target_time){
        done=true;
        dt=target_time-time;
        }*/
    T dt=simulation.dt;
    done=time+dt>=target_time;
    return dt;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void EVOLUTION<TV>::
Advance_One_Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    for(auto step: (*this)){step->Full_Step(simulation,dt,time);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void EVOLUTION<TV>::
Finalize()
{
    Json::Value root;
    for(auto step: (*this)){step->Finalize(root);}
    std::cout<<root<<std::endl;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EVOLUTION)
