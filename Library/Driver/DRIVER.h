//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class DRIVER
//##################################################################### 
#include <Evolution/EVOLUTION_COLLECTION.h>
#include <Driver/DRIVER.h>
#include <Data/DATA.h>
using namespace Mechanics;
//#####################################################################
template<class TV> DRIVER<TV>::
DRIVER(DATA<TV>& data)
    :data(data)
{
}
//#####################################################################
template<class TV> DRIVER<TV>::
~DRIVER()
{
}
//#####################################################################
template<class TV> void DRIVER<TV>::
Initialize()
{
    // for restarts, need to load data to current frame
    if(data.restart){
        data.Read_Output_Files(data.restart_frame);
        data.current_frame=data.restart_frame;
        data.Read_Time(data.current_frame);}
    output_number=data.current_frame;
    data.evolution->Update_Position_Based_State(data.time,*data.data); // make sure we're set for the current data
}
//#####################################################################
template<class TV> void DRIVER<TV>::
Advance_One_Time_Step(const T target_time,bool& done)
{
    // TODO: this became non-general.  Fix.  Generalize frame-not-frame concept.
    LOG_DEBUG::SCOPE scope("FRAME","Frame %d",data.current_frame+1);
    data.dt=data.evolution->Compute_Dt(data.exact_dt,data.minimum_dt,data.time,target_time,*data.data,done);
    LOG_DEBUG::cout<<"dt: "<<data.dt<<std::endl;
    data.evolution->Position_Step(data.dt,data.time,*data.data);
    data.evolution->Velocity_Step(data.dt,data.time,*data.data);
    data.current_frame++;
    LOG_DEBUG::cout<<"TIME = "<<data.time<<std::endl;
    data.time+=data.dt;
}
//#####################################################################
template<class TV> void DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    data.Write_Output_Files(data.current_frame);
    bool done=false;
    for(int substep=1;!done;substep++){
        Advance_One_Time_Step(target_time,done);
        data.Write_Output_Files(data.current_frame);
    }
}
//#####################################################################
template<class TV> void DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    if(!data.restart){
        data.Write_Output_Files(data.first_frame);
        data.Write_First_Frame(data.first_frame);}
    Advance_To_Target_Time(data.last_time);
}
//#####################################################################
template class DRIVER<VECTOR<float,1> >;
template class DRIVER<VECTOR<float,2> >;
template class DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DRIVER<VECTOR<double,1> >;
template class DRIVER<VECTOR<double,2> >;
template class DRIVER<VECTOR<double,3> >;
#endif
