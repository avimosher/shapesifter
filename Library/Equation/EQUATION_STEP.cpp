#include <Data/DATA.h>
#include <Equation/EQUATION.h>
#include <Equation/EQUATION_STEP.h>
#include <Evolution/QUALITY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> EQUATION_STEP<TV>::
EQUATION_STEP()
{}
///////////////////////////////////////////////////////////////////////
template<class TV> EQUATION_STEP<TV>::
~EQUATION_STEP()
{}
///////////////////////////////////////////////////////////////////////
template<class TV> void EQUATION_STEP<TV>::
Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    Matrix<T,Dynamic,1> positions;
    Matrix<T,Dynamic,1> velocities;
    //data.Pack_Velocities(velocities);velocities.setZero();
    data.Pack_Positions(positions);
    equation->Linearize(data,force,velocities,dt,time,true);
    //int count=0;
    Matrix<T,Dynamic,1> fractional_velocities;
    QUALITY solve_quality;
    while(!equation->Satisfied(data,force,fractional_velocities,dt,time)){
        // solve for variables
        fractional_velocities=equation->Solve(data,force,dt,time);
        // step data according to result
        data.Unpack_Positions(positions);
        if(!velocities.rows()){
            velocities.resize(fractional_velocities.rows(),1);
            velocities.setZero();
        }
        velocities+=(T).25*fractional_velocities;
        data.Unpack_Velocities(velocities.block(0,0,data.Velocity_DOF(),1));
        data.Unpack_Positions(positions);
        data.Step(solve_quality,fractional_velocities);

        /*Matrix<T,Dynamic,1> print_positions;
        data.Pack_Positions(print_positions);
        std::cout<<"Positions: "<<std::endl<<print_positions<<std::endl;
        std::cout<<"Velocities: "<<std::endl<<velocities<<std::endl;*/

        equation->Linearize(data,force,velocities,dt,time,false); // make force balance a force as well?
        
        //count++;if(count>80){exit(0);}
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EQUATION_STEP)
