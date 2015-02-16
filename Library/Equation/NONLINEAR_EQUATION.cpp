#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/EQUATION_STEP.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Evolution/EVOLUTION.h>
#include <Evolution/QUALITY.h>
#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <iostream>
#include <vector>
#include <Eigen/IterativeSolvers>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const Matrix<T,Dynamic,1>& velocities,const T dt,const T time,const bool stochastic)
{
    int full_size=data.size()+force.size();
    full_matrix.resize(full_size,full_size);
    full_right_hand_side.resize(full_size,1);

    std::vector<std::vector<Triplet<T>>> force_terms;
    full_right_hand_side(0,0).resize(data.Velocity_DOF(),1);
    full_right_hand_side(0,0).setZero();
    int size=full_right_hand_side(0,0).rows();
    J.resize(size,size);
    full_matrix(0,0).resize(size,size);
    // TODO: one block matrix per data type too
    force_terms.resize(data.size());
    for(int i=0;i<data.size();i++){
        data[i]->Inertia(force_terms[i]);
    }
    //for(int i=0;i<matrix.rows();i++){force_terms.push_back(Triplet<T>(i,i,1));} // identity portion
    for(int i=0;i<force.size();i++){
        // TODO: force_terms needs to be properly handled when there are multiple data types
        force[i]->Linearize(data,dt,time,force_terms[0],full_matrix(i+1,0),full_right_hand_side(0,0),full_right_hand_side(i+1,0),stochastic);
        //force[i]->Special(data,dt,time);
        full_matrix(0,i+1)=full_matrix(i+1,0).transpose();
    }
    // for the sake of sanity, assume that each force adds a constraint block as well as a contribution to the velocity block
    // Each such block will be required to be in terms of elementary T types, but they will remain separate.  This is a good compromise.
    // build matrix from force terms and constraint terms.  Not that this is sufficiently general...
    for(int i=0;i<data.size();i++){
        full_matrix(i,i).setFromTriplets(force_terms[i].begin(),force_terms[i].end());
    }
    Merge_Block_Matrices(full_matrix,J);
    //full_right_hand_side(0,0)-=full_matrix(0,0)*velocities;
    Merge_Block_Vectors(full_right_hand_side,right_hand_side);
    //right_hand_side_full=J*right_hand_side;
    //matrix=J.adjoint()*J;
    right_hand_side_full=right_hand_side;
    matrix=J;
    std::cout<<matrix<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar NONLINEAR_EQUATION<TV>::
Calculate_RHS_And_Norm(const DATA<TV>& data,const FORCE<TV>& force,const Matrix<T,Dynamic,1>& velocities)
{
    Matrix<T,Dynamic,1> forces;force.Pack_Forces(forces); // ask for stored forces appropriate for this solve
    Matrix<T,Dynamic,1> current_solution(right_hand_side.size());
    current_solution<<velocities,
        forces;
    int velocity_count=data.Velocity_DOF();
    right_hand_side.block(0,0,velocity_count,1)-=J.block(0,0,velocity_count,matrix.cols())*current_solution;
    //right_hand_side_full=J*right_hand_side;
    right_hand_side_full=right_hand_side;
    std::cout<<"Current solution: "<<current_solution.transpose()<<std::endl;
    std::cout<<"RHS: "<<right_hand_side.transpose()<<std::endl;
    return right_hand_side.squaredNorm();
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> NONLINEAR_EQUATION<TV>::
Solve(const Matrix<T,Dynamic,1>& guess)
{
    const int solve_iterations=200;
    MINRES<SparseMatrix<T>> solver;
    solver.setMaxIterations(solve_iterations);
    std::cout<<"Solve RHS: "<<right_hand_side_full.transpose()<<std::endl;
    solver.compute(matrix);
    auto solution=solver.solveWithGuess(right_hand_side_full,guess);
    std::cout<<"Solution: "<<solution.transpose()<<std::endl;
    return solution;
}
///////////////////////////////////////////////////////////////////////
template<class TV> bool NONLINEAR_EQUATION<TV>::
Satisfied(DATA<TV>& data,FORCE<TV>& force,const Matrix<T,Dynamic,1>& solve_result,QUALITY<T>& solve_quality)
{
    int velocity_count=data.Velocity_DOF();
    //auto residual=right_hand_side-matrix.block(0,velocity_count,matrix.rows(),matrix.cols()-velocity_count)*solve_result.block(velocity_count,0,solve_result.rows()-velocity_count,1);
    auto residual=right_hand_side-matrix*solve_result;
    T norm=residual.norm();
    solve_quality.Update(norm);
    //std::cout<<"solve result: "<<solve_result.transpose()<<std::endl;
    //std::cout<<"rhs: "<<right_hand_side.transpose()<<std::endl;
    //std::cout<<"residual: "<<residual.transpose()<<std::endl;
    std::cout<<"Norm: "<<norm<<std::endl;
    return norm<1e-6;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(NONLINEAR_EQUATION)
DEFINE_AND_REGISTER_PARSER(NONLINEAR_EQUATION,void)
{
    auto evolution_step=std::make_shared<EQUATION_STEP<TV>>();
    evolution_step->equation=new NONLINEAR_EQUATION<TV>();
    simulation.evolution.push_back(evolution_step);
    return 0;
}
