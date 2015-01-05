///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class EVOLUTION
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Equation/EQUATION.h>
#include <Evolution/EVOLUTION.h>
#include <Evolution/EVOLUTION_STEP.h>
#include <Evolution/EVOLUTION_TYPE.h>
#include <Evolution/QUALITY.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> EVOLUTION<TV>::
EVOLUTION()
{
    
}
///////////////////////////////////////////////////////////////////////
template<class TV> EVOLUTION<TV>::
~EVOLUTION()
{
    
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar EVOLUTION<TV>::
Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T time,const T target_time,bool& done)
{
    T dt=1;
    done=time+dt>=target_time;
    return dt;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void EVOLUTION<TV>::
Evolution_Step(EVOLUTION_STEP<TV>& step,DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    for(EVOLUTION_STEP<TV>* prerequisite_step : step.prerequisites) {
        if(!prerequisite_step->Satisfied(data,force,dt,time)) {
            Evolution_Step(*prerequisite_step,data,force,dt,time);
        }
    }
    EQUATION<TV>& equation=*step.equation;
    Matrix<T,Dynamic,1> positions;
    Matrix<T,Dynamic,1> velocities;
    //data.Pack_Velocities(velocities);velocities.setZero();
    data.Pack_Positions(positions);
    equation.Linearize(data,force,velocities,dt,time,true);
    //int count=0;
    while(!equation.Satisfied(data,force,dt,time)){
        // solve for variables
        Matrix<T,Dynamic,1> fractional_velocities=equation.Solve(data,force,dt,time);
        QUALITY solve_quality;
        // step data according to result
        data.Unpack_Positions(positions);
        if(!velocities.rows()){
            velocities.resize(fractional_velocities.rows(),1);
            velocities.setZero();
        }
        velocities+=(T).5*fractional_velocities;
        data.Unpack_Velocities(velocities.block(0,0,data.Velocity_DOF(),1));
        data.Unpack_Positions(positions);
        data.Step(solve_quality,fractional_velocities);

        /*Matrix<T,Dynamic,1> print_positions;
        data.Pack_Positions(print_positions);
        std::cout<<"Positions: "<<std::endl<<print_positions<<std::endl;
        std::cout<<"Velocities: "<<std::endl<<velocities<<std::endl;*/

        equation.Linearize(data,force,velocities,dt,time,false); // make force balance a force as well?
        
        //count++;if(count>20){exit(0);}
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void EVOLUTION<TV>::
Advance_One_Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    for(auto step : (*this)) {
        if(!step->Satisfied(data,force,dt,time)) {
            Evolution_Step(*step,data,force,dt,time);
        }
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EVOLUTION)
