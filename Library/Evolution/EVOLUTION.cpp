//---------------------------------------------------------------------
// Copyright 2014, Avi Robinson-Mosher.
//---------------------------------------------------------------------
// Class EVOLUTION
//---------------------------------------------------------------------
#include <Evolution/EVOLUTION.h>
using namespace Mechanics;
//---------------------------------------------------------------------
template<class TV> EVOLUTION<TV>::
EVOLUTION()
{
}
//---------------------------------------------------------------------
template<class TV> typename TV::Scalar EVOLUTION<TV>::
Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time)
{

}
//---------------------------------------------------------------------
template<class TV> void EVOLUTION<TV>::
Position_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    for(EVOLUTION_TYPE<TV>& evolution : evolution_collection) {
        evolution.Position_Step(data,force,dt,time);
    }
    
}
//---------------------------------------------------------------------
template<class TV> void EVOLUTION<TV>::
Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    force.Compute(data,dt,time);
    Explicit_Velocity_Step(data,force,dt,time);
    Implicit_Velocity_Step(data,force,dt,time);
}
GENERIC_TYPE_DEFINITION(EVOLUTION)
