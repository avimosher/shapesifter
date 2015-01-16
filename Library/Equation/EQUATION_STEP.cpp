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
    Matrix<T,Dynamic,1> current_velocities;
    Matrix<T,Dynamic,1> current_forces;
    Matrix<T,Dynamic,1> solve_forces,solve_velocities; // TODO: make these blocks?
    current_velocities.resize(data.Velocity_DOF(),1);current_velocities.setZero();
    data.Pack_Positions(positions);
    equation->Linearize(data,force,current_velocities,dt,time,true);

    //data.Pack_Velocities(current_solution);
    //current_solution.setZero();
    force.Pack_Forces(solve_forces);
    solve_velocities=current_velocities;
    Matrix<T,Dynamic,1> solve_vector;solve_vector.resize(solve_velocities.size()+solve_forces.size());
    solve_vector<<solve_velocities,
        solve_forces;

    //int count=0;
    QUALITY solve_quality;
    while(!equation->Satisfied(data,force,solve_vector,dt,time)){
        // solve for variables
        solve_vector=equation->Solve(solve_vector);
        solve_velocities=solve_vector.block(0,0,data.Velocity_DOF(),1);
        solve_forces=solve_vector.block(data.Velocity_DOF(),0,solve_vector.rows()-data.Velocity_DOF(),1);

        current_velocities+=(T).25*solve_velocities;

        // NOTE: for the sake of things like snap constraints, it's good that this puts the velocity in data (and the initial velocity should probably be zero)
        // step data according to result
        data.Unpack_Positions(positions);
        data.Unpack_Velocities(current_velocities);
        data.Unpack_Positions(positions);
        data.Step(solve_quality,current_velocities);

        /*Matrix<T,Dynamic,1> print_positions;
        data.Pack_Positions(print_positions);
        std::cout<<"Positions: "<<std::endl<<print_positions<<std::endl;
        std::cout<<"Velocities: "<<std::endl<<velocities<<std::endl;*/

        force.Unpack_Forces(solve_forces);
        equation->Linearize(data,force,current_velocities,dt,time,false); // make force balance a force as well?
        solve_velocities.setZero();
        force.Pack_Forces(solve_forces);
        solve_vector.resize(solve_velocities.size()+solve_forces.size());
        solve_vector<<solve_velocities,
            solve_forces;
        //count++;if(count>80){exit(0);}
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EQUATION_STEP)
