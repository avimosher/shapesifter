///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class EVOLUTION
///////////////////////////////////////////////////////////////////////
#include <Evolution/EVOLUTION.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> EVOLUTION<TV>::
EVOLUTION()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar EVOLUTION<TV>::
Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time)
{

}
///////////////////////////////////////////////////////////////////////
template<class TV> void EVOLUTION<TV>::
Position_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    // for inertial:
    // advance position with current-ish velocities (possibly update to t+1/2 first)
    // so, compute forces to correct time, do velocity solve (enforce constraints)
    // for non-inertial:
    // compute forces to correct time, do position solve (enforce constraints)
    // solve must operate on correct unknowns

    EQUATION<TV> equation;
    do{
        // sets up structure but does not solve for variables (forces and data)
        equation.Linearize(data,force,dt,time); // make force balance a force as well?
        // solve for variables
        equation.Solve(data,force,dt,time);
    } while(!equation.Satisfied(data,force,dt,time));
}
///////////////////////////////////////////////////////////////////////
template<class TV> void EVOLUTION<TV>::
Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    force.Compute(data,dt,time);
    Explicit_Velocity_Step(data,force,dt,time);
    Implicit_Velocity_Step(data,force,dt,time);
}
GENERIC_TYPE_DEFINITION(EVOLUTION)
