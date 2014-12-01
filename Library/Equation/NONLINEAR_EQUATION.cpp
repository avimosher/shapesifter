///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class NONLINEAR_EQUATION
///////////////////////////////////////////////////////////////////////
#include <Evolution/NONLINEAR_EQUATION.h>
#include <Eigen/IterativeSolvers>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> NONLINEAR_EQUATION<TV>::
NONLINEAR_EQUATION()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    Triplet<T> force_terms;
    Matrix<T,Dynamic,1> right_hand_side;
    for(FORCE_TYPE<TV>* force_type : force){
        // Eigen nicely sums duplicate entries in a Triplet list - perfect.
        Triple<T> constraint_terms;
        force_type.Linearize(data,dt,time,force_terms,constraint_terms,right_hand_side);
    }
    // build matrix fromo force terms and constraint terms.  Not that this is sufficiently general...
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> NONLINEAR_EQUATION<TV>::
Solve(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    const int newton_iterations=40;
    const int solve_iterations=200;
    MINRES<SparseMatrix<T> > solver;
    solver.setMaxIterations(solve_iterations);

    for(int i=0;i<newton_iterations;i++){
        solver.compute(A);
        x=solver.solve(b);
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(NONLINEAR_EQUATION)
