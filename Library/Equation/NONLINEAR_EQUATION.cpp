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
template<class TV> typename TV::Scalar NONLINEAR_EQUATION<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic)
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
        data[i]->Inertia(force_terms[i],full_right_hand_side[i]);
    }
    for(int i=0;i<force.size();i++){
        // TODO: force_terms needs to be properly handled when there are multiple data types
        force[i]->Linearize(data,dt,time,force_terms[0],full_matrix(i+1,0),full_right_hand_side(0,0),full_right_hand_side(i+1,0),stochastic);
        full_matrix(0,i+1)=full_matrix(i+1,0).transpose();
    }
    // for the sake of sanity, assume that each force adds a constraint block as well as a contribution to the velocity block
    // Each such block will be required to be in terms of elementary T types, but they will remain separate.  This is a good compromise.
    // build matrix from force terms and constraint terms.  Not that this is sufficiently general...
    for(int i=0;i<data.size();i++){
        full_matrix(i,i).setFromTriplets(force_terms[i].begin(),force_terms[i].end());
    }
    Merge_Block_Matrices(full_matrix,J);
    Merge_Block_Vectors(full_right_hand_side,right_hand_side);
    //right_hand_side_full=J*right_hand_side;
    //matrix=J.adjoint()*J;
    right_hand_side_full=right_hand_side;
    matrix=J;
    Matrix<T,Dynamic,1> rowExtrema(matrix.rows());rowExtrema.setZero();
    for(int k=0;k<matrix.outerSize();k++){
        for(typename SparseMatrix<T>::InnerIterator it(matrix,k);it;++it){
            rowExtrema[it.row()]=std::max(rowExtrema[it.row()],std::abs(it.value()));
        }
    }
    conditioner=rowExtrema;//.cwiseSqrt().cwiseInverse();
    //std::cout<<"Conditioner: "<<conditioner.transpose()<<std::endl;
    //matrix=conditioner.asDiagonal()*matrix*conditioner.asDiagonal();
    //right_hand_side_full=right_hand_side_full;
    std::cout<<matrix<<std::endl;
    std::cout<<"RHS: "<<right_hand_side_full.transpose()<<std::endl;
    return right_hand_side_full.squaredNorm();
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> NONLINEAR_EQUATION<TV>::
Solve()
{
    const int solve_iterations=200;
    //MINRES<SparseMatrix<T>,Lower,RowPreconditioner<T>> solver;
    MINRES<SparseMatrix<T>> solver;
    solver.compute(matrix);
    //solver.preconditioner().SetDiagonal(conditioner);
    solver.setMaxIterations(solve_iterations);
    Matrix<T,Dynamic,1> solution=solver.solve(right_hand_side_full);
    //std::cout<<"Iterations: "<<solver.iterations()<<std::endl;
    std::cout<<"Solution: "<<solution.transpose()<<std::endl;
    std::cout<<"Solve RHS: "<<right_hand_side_full.transpose()<<std::endl;
    std::cout<<"A*x: "<<(matrix*solution).transpose()<<std::endl;
    return solution;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar NONLINEAR_EQUATION<TV>::
Sufficient_Descent_Factor(const Matrix<T,Dynamic,1>& direction)
{
    return -direction.dot(matrix*right_hand_side_full);
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
