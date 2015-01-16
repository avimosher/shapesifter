#include <Data/DATA.h>
#include <Equation/EQUATION.h>
#include <Equation/EQUATION_STEP.h>
#include <Evolution/QUALITY.h>
#include <Force/FORCE.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void EQUATION_STEP<TV>::
Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    Matrix<T,Dynamic,1> positions;
    Matrix<T,Dynamic,1> current_solution;
    //data.Pack_Velocities(velocities);velocities.setZero();
    data.Pack_Positions(positions);
    equation->Linearize(data,force,current_solution,dt,time,true);

    data.Pack_Velocities(current_solution);
    current_solution.setZero();
    force.Pack_Forces(current_solution);

    Matrix<T,Dynamic,1> solve_vector;

    //int count=0;
    QUALITY solve_quality;
    while(!equation->Satisfied(data,force,solve_vector,dt,time)){
        // solve for variables
        solve_vector=equation->Solve(solve_vector);
        // step data according to result
        data.Unpack_Positions(positions);
        current_solution+=(T).25*solve_vector;
        // NOTE: for the sake of things like snap constraints, it's good that this puts the velocity in data (and the initial velocity should probably be zero)
        data.Unpack_Velocities(current_solution.block(0,0,data.Velocity_DOF(),1));
        data.Unpack_Positions(positions);
        data.Step(solve_quality,current_solution);

        /*Matrix<T,Dynamic,1> print_positions;
        data.Pack_Positions(print_positions);
        std::cout<<"Positions: "<<std::endl<<print_positions<<std::endl;
        std::cout<<"Velocities: "<<std::endl<<velocities<<std::endl;*/

        force.Unpack_Forces(solve_vector);
        equation->Linearize(data,force,current_solution,dt,time,false); // make force balance a force as well?
        solve_vector.setZero();
        force.Pack_Forces(solve_vector);
        //count++;if(count>80){exit(0);}
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EQUATION_STEP)
