#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/EQUATION.h>
#include <Equation/EQUATION_STEP.h>
#include <Evolution/QUALITY.h>
#include <Force/FORCE.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void EQUATION_STEP<TV>::
Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    Matrix<T,Dynamic,1> positions;
    Matrix<T,Dynamic,1> current_velocities;
    Matrix<T,Dynamic,1> current_forces;
    Matrix<T,Dynamic,1> solve_forces,solve_velocities; // TODO: make these blocks?
    DATA<TV>& data=simulation.data;
    FORCE<TV>& force=simulation.force;
    current_velocities.resize(data.Velocity_DOF(),1);current_velocities.setZero();
    data.Pack_Positions(positions);
    equation->Linearize(data,force,current_velocities,dt,time,true);

    force.Pack_Forces(solve_forces);
    solve_velocities=current_velocities;
    Matrix<T,Dynamic,1> solve_vector;solve_vector.resize(solve_velocities.size()+solve_forces.size());
    solve_vector<<solve_velocities,
        solve_forces;

    int count=0;
    QUALITY solve_quality;
    while(!equation->Satisfied(data,force,solve_vector,dt,time)){
        solve_vector=equation->Solve(solve_vector);
        solve_velocities=solve_vector.block(0,0,data.Velocity_DOF(),1);
        solve_forces=solve_vector.block(data.Velocity_DOF(),0,solve_vector.rows()-data.Velocity_DOF(),1);

        current_velocities+=(T).15*solve_velocities;

        // NOTE: for the sake of things like snap constraints, it's good that this puts the velocity in data (and the initial velocity should probably be zero)
        // step data according to result
        data.Unpack_Positions(positions);
        data.Unpack_Velocities(current_velocities);
        data.Unpack_Positions(positions);
        data.Step(solve_quality);

        force.Unpack_Forces(solve_forces);
        equation->Linearize(data,force,current_velocities,dt,time,false); // make force balance a force as well?
        solve_velocities.setZero();
        force.Pack_Forces(solve_forces);
        solve_vector.resize(solve_velocities.size()+solve_forces.size());
        solve_vector<<solve_velocities,
            solve_forces;
        count++;
        if(simulation.substeps){
            simulation.title="Substep for "+std::to_string(count)+" frame "+std::to_string(simulation.current_frame);
            simulation.Write(simulation.current_frame);
            simulation.current_frame++;
            std::cout<<"Write frame: "<<simulation.current_frame<<std::endl;
        }
        //count++;if(count>35){exit(0);}
    }
    std::cout<<"Steps: "<<count<<std::endl;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EQUATION_STEP)
