#include <Driver/SIMULATION.h>
#include <Equation/LINE_SEARCH.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> LINE_SEARCH<TV>::
LINE_SEARCH()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> void LINE_SEARCH<TV>::
Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
/*
    Matrix<T,Dynamic,1> positions;
    Matrix<T,Dynamic,1> current_velocities;
    Matrix<T,Dynamic,1> current_forces;
    Matrix<T,Dynamic,1> solve_velocities; // TODO: make these blocks?
    STORED_FORCE<T> solve_forces;
    DATA<TV>& data=simulation.data;
    FORCE<TV>& force=simulation.force;
    current_velocities.resize(data.Velocity_DOF(),1);current_velocities.setZero();
    data.Unpack_Velocities(current_velocities);
    data.Pack_Positions(positions);
    force.Pack_Forces(solve_forces);
    solve_forces.setZero();
    force.Unpack_Forces(solve_forces);
    equation->Linearize(data,force,dt,time,true);
    T last_norm=equation->Evaluate();

    std::cout<<"First norm: "<<last_norm<<std::endl;

    force.Pack_Forces(solve_forces);
    solve_velocities=current_velocities;

    int count=0;
    QUALITY<T> solve_quality;
    T epsilon=1e-8;
    int step_limit=16;
    T c1=.5,c2=.9;
    T last_ratio=.5;
    while(last_norm>epsilon){
        auto solve_vector=equation->Solve();
        solve_velocities=solve_vector.block(0,0,data.Velocity_DOF(),1);
        solve_forces.Set(solve_vector.block(data.Velocity_DOF(),0,solve_vector.rows()-data.Velocity_DOF(),1));

        T sufficient_descent_factor=equation->Sufficient_Descent_Factor(solve_vector);
        T ratio=2*last_ratio,norm;
        int i=0;

        for(;i<step_limit&&ratio>epsilon;i++){
            // update state
            data.Unpack_Positions(positions); 
            force.Increment_Forces(solve_forces,ratio);
            data.Unpack_Velocities(current_velocities+ratio*solve_velocities);
            data.Step();
            if(simulation.substeps){
                simulation.Write("Frame "+std::to_string(simulation.current_frame)+" substep "+std::to_string(count)+" line search "+std::to_string(i));
                std::cout<<"Write frame: "<<simulation.current_frame<<std::endl;
            }
            equation->Linearize(data,force,dt,time,false); // linearize around a point and calculate norm there
            norm=equation->Evaluate();
            std::cout<<"Norm with ratio "<<ratio<<" is "<<norm<<" ("<<last_norm<<")"<<std::endl;
            //T curvature_factor=equation->Sufficient_Descent_Factor(solve_vector);
            if(norm<=last_norm+c1*sufficient_descent_factor*ratio){// && curvature_factor<0 && std::abs(curvature_factor)<=std::abs(c2*sufficient_descent_factor)){ //TODO: properly implement "sufficient reduction" criterion
                break;
            }
            // restore previous state
            force.Increment_Forces(solve_forces,-ratio);
            ratio/=2;
        }
        if(i==step_limit || ratio<=epsilon){std::cout<<"STALLED"<<std::endl;break;}
        last_ratio=ratio;
        current_velocities+=ratio*solve_velocities;
        last_norm=norm;
        force.Pack_Forces(solve_forces); // this resizes the forces correctly

        count++;
    }
    std::cout<<"Steps: "<<count<<std::endl;
*/
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(LINE_SEARCH)
DEFINE_AND_REGISTER_PARSER(LINE_SEARCH,void)
{
    auto step=std::make_shared<LINE_SEARCH<TV>>();
    step->equation=std::make_shared<NONLINEAR_EQUATION<TV>>();
    simulation.evolution.push_back(step);
    return 0;
}

