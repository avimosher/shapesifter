///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class DRIVER
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Driver/DRIVER.h>
#include <Evolution/EVOLUTION.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> DRIVER<TV>::
DRIVER(DATA<TV>& data,EVOLUTION<TV>& evolution)
    :data(data),evolution(evolution)
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> DRIVER<TV>::
~DRIVER()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Initialize()
{
    // for restarts, need to load data to current frame
    if(data.restart){
        data.Read_Output_Files(restart_frame);
        current_frame=restart_frame;
        data.Read_Time(current_frame);}
    output_number=current_frame;
    evolution.Update_Position_Based_State(time,data); // make sure we're set for the current data
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Advance_One_Time_Step(const T target_time,bool& done)
{
    // TODO: this became non-general.  Fix.  Generalize frame-not-frame concept.
    LOG_DEBUG::SCOPE scope("FRAME","Frame %d",data.current_frame+1);
    T dt=evolution.Compute_Dt(data,force,target_time);
    LOG_DEBUG::cout<<"dt: "<<dt<<std::endl;
    evolution.Position_Step(data,force,dt,time);
    evolution.Velocity_Step(data,force,dt,time);
    current_frame++;
    LOG_DEBUG::cout<<"TIME = "<<time<<std::endl;
    time+=dt;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    data.Write_Output_Files(current_frame);
    bool done=false;
    for(int substep=1;!done;substep++){
        Advance_One_Time_Step(target_time,done);
        data.Write_Output_Files(current_frame);
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    if(!restart){
        data.Write_Output_Files(first_frame);
        data.Write_First_Frame(first_frame);}
    Advance_To_Target_Time(last_time);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(DRIVER)
