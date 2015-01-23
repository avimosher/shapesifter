#include <Data/DATA.h>
#include <Driver/DRIVER.h>
#include <Driver/SIMULATION.h>
#include <Evolution/EVOLUTION.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Initialize()
{
    if(simulation->restart){
        simulation->current_frame=simulation->restart_frame;
        simulation->Read(simulation->restart_frame);
    }
    //evolution.Update_Position_Based_State(time,data); // make sure we're set for the current data
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Advance_One_Time_Step(const T target_time,bool& done)
{
    // TODO: this became non-general.  Fix.  Generalize frame-not-frame concept.
    T dt=simulation->evolution.Compute_Dt(simulation->data,simulation->force,simulation->time,target_time,done);
    simulation->evolution.Advance_One_Step(simulation->data,simulation->force,dt,simulation->time);
    simulation->current_frame++;
    simulation->time+=dt;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        Advance_One_Time_Step(target_time,done);
        simulation->Write(simulation->current_frame);
        std::cout<<simulation->current_frame<<std::endl;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Advance_To_Target_Time(simulation->last_time);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(DRIVER)

