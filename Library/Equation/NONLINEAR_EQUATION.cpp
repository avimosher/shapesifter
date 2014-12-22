///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class NONLINEAR_EQUATION
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Evolution/EVOLUTION.h>
#include <Evolution/EVOLUTION_STEP.h>
#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <iostream>
#include <vector>
#include <Eigen/IterativeSolvers>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    std::vector<Triplet<T>> force_terms;
    data.Variables(right_hand_side);
    int size=right_hand_side.rows();
    matrix.resize(size,size);
    for(int i=0;i<matrix.rows();i++){
        force_terms.push_back(Triplet<T>(i,i,1));
    }
    for(auto& force_type : force){
        // Eigen nicely sums duplicate entries in a Triplet list - perfect.
        std::vector<Triplet<T>> constraint_terms;
        force_type->Linearize(data,dt,time,force_terms,constraint_terms,right_hand_side);
    }
    // build matrix from force terms and constraint terms.  Not that this is sufficiently general...
    matrix.setFromTriplets(force_terms.begin(),force_terms.end());
    std::cout<<matrix<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> NONLINEAR_EQUATION<TV>::
Solve(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    const int solve_iterations=200;
    MINRES<SparseMatrix<T> > solver;
    solver.setMaxIterations(solve_iterations);
    solver.compute(matrix);
    return solver.solve(right_hand_side);
}
///////////////////////////////////////////////////////////////////////
template<class TV> bool NONLINEAR_EQUATION<TV>::
Satisfied(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    Matrix<T,Dynamic,1> x;
    data.Variables(x); // TODO: this should probably be per solve type
    return (matrix*x-right_hand_side).squaredNorm()<1e-8;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(NONLINEAR_EQUATION)
DEFINE_AND_REGISTER_PARSER(NONLINEAR_EQUATION)
{
    auto evolution_step=std::make_shared<EVOLUTION_STEP<TV>>();
    evolution_step->equation=new NONLINEAR_EQUATION<TV>();
    simulation.evolution.push_back(evolution_step);
}
