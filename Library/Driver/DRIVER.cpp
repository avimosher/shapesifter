#include <Data/DATA.h>
#include <Driver/DRIVER.h>
#include <Driver/SIMULATION.h>
#include <Evolution/EVOLUTION.h>
#include <Utilities/LOG.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Initialize()
{
    if(simulation->restart){
        simulation->current_frame=simulation->restart_frame;
        simulation->output_number=simulation->restart_frame;
        simulation->Read(simulation->restart_frame-1);
    }
    else{simulation->Write("Initial state");}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Advance_One_Time_Step(const T target_time,bool& done)
{
    // TODO: this became non-general.  Fix.  Generalize frame-not-frame concept.
    T dt=simulation->evolution.Compute_Dt(*simulation,simulation->time,target_time,done);
    simulation->evolution.Advance_One_Step(*simulation,dt,simulation->time);
    simulation->time+=dt;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        Advance_One_Time_Step(target_time,done);
        //simulation->Write("End frame "+std::to_string(simulation->current_frame));
        simulation->current_frame++;
        LOG::cout<<"\n\nFrame "<<simulation->current_frame<<std::endl;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Advance_To_Target_Time(simulation->last_time);
    simulation->evolution.Finalize();
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(DRIVER)

