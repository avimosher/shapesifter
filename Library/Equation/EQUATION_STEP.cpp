#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/EQUATION.h>
#include <Equation/EQUATION_STEP.h>
#include <Equation/NONLINEAR_EQUATION.h>
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
    data.Unpack_Velocities(current_velocities);
    data.Pack_Positions(positions);
    T last_norm=equation->Linearize(data,force,dt,time,true);
    force.Pack_Forces(solve_forces);
    solve_forces.setZero();
    force.Unpack_Forces(solve_forces);

    std::cout<<"First norm: "<<last_norm<<std::endl;

    force.Pack_Forces(solve_forces);
    solve_velocities=current_velocities;

    int count=0;
    QUALITY<T> solve_quality;
    T epsilon=1e-8;
    int step_limit=16;
    T c1=.5,c2=.9;
    while(last_norm>epsilon){
        auto solve_vector=equation->Solve();
        //std::cout<<"Solve vector: "<<solve_vector.transpose()<<std::endl;
        solve_velocities=solve_vector.block(0,0,data.Velocity_DOF(),1);
        solve_forces=solve_vector.block(data.Velocity_DOF(),0,solve_vector.rows()-data.Velocity_DOF(),1);

        T sufficient_descent_factor=equation->Sufficient_Descent_Factor(solve_vector);
        std::cout<<"Sufficient descent: "<<sufficient_descent_factor<<std::endl;
        T ratio=1,norm;
        int i=0;

        for(;i<step_limit;i++){
            // update state
            data.Unpack_Positions(positions); 
            force.Increment_Forces(ratio*solve_forces);
            data.Unpack_Velocities(current_velocities+ratio*solve_velocities);
            data.Step();
            norm=equation->Linearize(data,force,dt,time,false); // linearize around a point and calculate norm there
            std::cout<<"Norm with ratio "<<ratio<<" is "<<norm<<" ("<<last_norm<<")"<<std::endl;
            T curvature_factor=equation->Sufficient_Descent_Factor(solve_vector);
            std::cout<<"Curvature factor: "<<curvature_factor<<std::endl;
            if(norm<=last_norm+c1*sufficient_descent_factor*ratio && curvature_factor<0 && abs(curvature_factor)<=abs(c2*sufficient_descent_factor)){ //TODO: properly implement "sufficient reduction" criterion
                break;
            }
            // restore previous state
            force.Increment_Forces(-ratio*solve_forces);
            ratio/=2;
        }
        if(i==step_limit){break;}
        current_velocities+=ratio*solve_velocities;
        last_norm=norm;

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
