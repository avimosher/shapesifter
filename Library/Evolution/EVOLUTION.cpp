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
    for(auto step: (*this)){step->Finalize();}
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EVOLUTION)
