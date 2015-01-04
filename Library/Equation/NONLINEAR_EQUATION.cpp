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
#include <Utilities/EIGEN_HELPERS.h>
#include <iostream>
#include <vector>
#include <Eigen/IterativeSolvers>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const Matrix<T,Dynamic,1>& velocities,const T dt,const T time,const bool stochastic)
{
    int full_size=1+force.size();
    full_matrix.resize(full_size,full_size);
    full_right_hand_side.resize(full_size,1);

    std::vector<Triplet<T>> force_terms;
    full_right_hand_side(0,0).resize(data.Velocity_DOF(),1);
    full_right_hand_side(0,0).setZero();
    int size=full_right_hand_side(0,0).rows();
    matrix.resize(size,size);
    full_matrix(0,0).resize(size,size);
    for(int i=0;i<matrix.rows();i++){
        force_terms.push_back(Triplet<T>(i,i,1));
    }
    for(int i=0;i<force.size();i++){
        // Eigen nicely sums duplicate entries in a Triplet list - perfect.
        force[i]->Linearize(data,dt,time,force_terms,full_matrix(i+1,0),full_right_hand_side(0,0),full_right_hand_side(i+1,0),stochastic);
        full_matrix(0,i+1)=full_matrix(i+1,0).transpose();
    }
    //std::cout<<full_right_hand_side(0,0)<<std::endl;
    // for the sake of sanity, assume that each force adds a constraint block as well as a contribution to the velocity block
    // Each such block will be required to be in terms of elementary T types, but they will remain separate.  This is a good compromise.
    // build matrix from force terms and constraint terms.  Not that this is sufficiently general...
    full_matrix(0,0).setFromTriplets(force_terms.begin(),force_terms.end());
    //full_right_hand_side(0,0)-=full_matrix(0,0)*velocities;
    Merge_Block_Matrices(full_matrix,matrix);
    if(velocities.rows()){
        full_right_hand_side(0,0)-=matrix.block(0,0,data.Velocity_DOF(),matrix.cols())*velocities;
    }
    Merge_Block_Vectors(full_right_hand_side,right_hand_side);
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> NONLINEAR_EQUATION<TV>::
Solve(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    const int solve_iterations=200;
    MINRES<SparseMatrix<T>> solver;
    solver.setMaxIterations(solve_iterations);
    solver.compute(matrix);
    // TODO: return only the velocity part?  Probably not.
    return solver.solve(right_hand_side);//.block(0,0,full_matrix(0,0).size(),1);
    //return right_hand_side(0,0);
}
///////////////////////////////////////////////////////////////////////
template<class TV> bool NONLINEAR_EQUATION<TV>::
Satisfied(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    std::cout<<"Satisfaction: "<<right_hand_side.squaredNorm()<<std::endl;
    return right_hand_side.squaredNorm()<1e-8;
    /*Matrix<T,Dynamic,1> x;
    data.Variables(x); // TODO: this should probably be per solve type
    return (matrix*x-right_hand_side).squaredNorm()<1e-8;*/
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(NONLINEAR_EQUATION)
DEFINE_AND_REGISTER_PARSER(NONLINEAR_EQUATION)
{
    auto evolution_step=std::make_shared<EVOLUTION_STEP<TV>>();
    evolution_step->equation=new NONLINEAR_EQUATION<TV>();
    simulation.evolution.push_back(evolution_step);
}
