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
    data.Pack_Positions(positions);
    equation->Linearize(data,force,current_velocities,dt,time,true);
    T last_norm=equation->Calculate_RHS_And_Norm(data,force,current_velocities);
    std::cout<<"First norm: "<<last_norm<<std::endl;

    force.Pack_Forces(solve_forces);
    solve_velocities=current_velocities;
    Matrix<T,Dynamic,1> solve_vector;solve_vector.resize(solve_velocities.size()+solve_forces.size());
    solve_vector<<solve_velocities,
        solve_forces;

    int count=0;
    QUALITY<T> solve_quality;
    T lambda=5;
    T epsilon=1e-8;
    int bad_count=0;
    while(last_norm>epsilon){//!equation->Satisfied(data,force,solve_vector,solve_quality)){
        std::cout<<"Entering loop"<<std::endl;
        solve_vector=equation->Solve(solve_vector,lambda);
        solve_velocities=solve_vector.block(0,0,data.Velocity_DOF(),1);

        // setting up positions; forces handle themselves
        data.Unpack_Positions(positions);
        data.Unpack_Velocities(current_velocities+solve_velocities);
        data.Step();
        solve_forces=solve_vector.block(data.Velocity_DOF(),0,solve_vector.rows()-data.Velocity_DOF(),1);
        force.Increment_Forces(solve_forces);
        std::cout<<"Before first loop linearize"<<std::endl;
        equation->Linearize(data,force,current_velocities+solve_velocities,dt,time,false); // linearize around a point and calculate norm there
        T norm=equation->Calculate_RHS_And_Norm(data,force,current_velocities+solve_velocities);

        T ratio=1;
        std::cout<<"Norm with ratio "<<ratio<<" is "<<norm<<", last norm "<<last_norm<<std::endl;

        while(norm>last_norm){
            force.Increment_Forces(-ratio*solve_forces);
            ratio/=2;
            data.Unpack_Positions(positions);
            data.Unpack_Velocities(current_velocities+ratio*solve_velocities);
            data.Step();
            force.Increment_Forces(ratio*solve_forces);
            equation->Linearize(data,force,current_velocities+solve_velocities,dt,time,false); // linearize around a point and calculate norm there
            norm=equation->Calculate_RHS_And_Norm(data,force,current_velocities+ratio*solve_velocities);
            std::cout<<"Norm with ratio "<<ratio<<" is "<<norm<<std::endl;
        }
        current_velocities+=ratio*solve_velocities;
        last_norm=norm;


#if 0
        // trust region
        solve_vector=equation->Solve(solve_vector,lambda);
        solve_velocities=solve_vector.block(0,0,data.Velocity_DOF(),1);

        // setting up positions; forces handle themselves
        data.Unpack_Positions(positions);
        data.Unpack_Velocities(current_velocities+solve_velocities);
        data.Step();

        solve_forces=solve_vector.block(data.Velocity_DOF(),0,solve_vector.rows()-data.Velocity_DOF(),1);
        force.Increment_Forces(solve_forces);

        // agh, I think this requires delta force too...
        std::cout<<"LINEARIZE AROUND NEW POSITIONS"<<std::endl;
        equation->Linearize(data,force,current_velocities+solve_velocities,dt,time,false); // linearize around a point and calculate norm there
        T norm=equation->Calculate_RHS_And_Norm(data,force,current_velocities+solve_velocities);
        T rho=(last_norm-norm)/(last_norm); // unless I'm greatly mistaken, my estimate for m(nu)-m(0) will be last_norm, since I'm solving for...no, that's not entirely true; it depends
                                            // on the tolerance on my solve.  Hmm.  Pretend it's like this for the moment and adjust later.

        std::cout<<"Solve vector: "<<solve_vector.transpose()<<std::endl;
        std::cout<<"Rho: "<<rho<<" lambda: "<<lambda<<" current norm: "<<norm<<" last norm: "<<last_norm<<std::endl;
        if(rho<0){
            // reject this step and reduce rho
            lambda*=4;
            data.Unpack_Positions(positions);
            std::cout<<"LINEARIZE TO FIX"<<std::endl;
            equation->Linearize(data,force,current_velocities,dt,time,false);
            force.Increment_Forces(-solve_forces); // this will work because the force set will be right...except I need to get this properly into the RHS
            equation->Calculate_RHS_And_Norm(data,force,current_velocities);
        }
        else{
            if(rho>.75){
                // doin' pretty good.  Increase lambda.
                lambda/=2;
            }
            //else if(rho<.25){
            //    lambda/=4;
            //}
            current_velocities+=solve_velocities;
            last_norm=norm;
        }
#endif
        /*

        T m=solve_vector.transpose()*equation->Gradient();
        T alpha;
        T initial_norm;// TODO must be set
        T t=-c*m;
        solve_velocities=solve_vector.block(0,0,data.Velocity_DOF(),1);


        // line search
        while(initial_norm-equation->Norm(alpha*solve_vector)>alpha*t){
            // update positions properly for linearization
            data.Unpack_Positions(positions);
            data.Unpack_Velocities(current_velocities+alpha*solve_velocities); // get positions from these velocities
            data.Step();
            force.Unpack_Forces(
            equation->Linearize(data,force,alpha*solve_velocities,dt,time,false); // linearize equations at candidate positions
            // assume new velocities/forces are zero and calculate norm (which is based on rhs)

            // if converged, can leave the result alone?  But need to increment forces.  Also, will now be linearized around correct point for next step!
            }*/


        /*current_velocities+=alpha*solve_velocities;
        force.Increment_Forces(solve_vector.block(data.Velocity_DOF(),0,solve_vector.rows()-data.Velocity_DOF(),1));
        //solve_forces=solve_vector.block(data.Velocity_DOF(),0,solve_vector.rows()-data.Velocity_DOF(),1);
        
        std::cout<<"Solve scale: "<<solve_quality.Scale_Result()<<std::endl;
        //current_velocities+=solve_quality.Scale_Result()*solve_velocities;

        // NOTE: for the sake of things like snap constraints, it's good that this puts the velocity in data (and the initial velocity should probably be zero)
        // step data according to result
        //data.Unpack_Positions(positions);
        data.Unpack_Velocities(current_velocities);
        data.Unpack_Positions(positions);
        data.Step();

        force.Unpack_Forces(solve_forces);
        equation->Linearize(data,force,current_velocities,dt,time,false); // make force balance a force as well?
        solve_velocities.setZero();
        force.Pack_Forces(solve_forces);
        solve_vector.resize(current_velocities.size()+solve_forces.size());solve_vector.setZero();
        //solve_vector<<solve_velocities,
        //    solve_forces;

        */
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
