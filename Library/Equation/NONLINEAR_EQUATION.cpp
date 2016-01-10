#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/MATRIX_BUNDLE.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Evolution/EVOLUTION.h>
#include <Evolution/QUALITY.h>
#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/LOG.h>
#include <iostream>
#include <vector>
#include <Eigen/IterativeSolvers>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> NONLINEAR_EQUATION<TV>::
Get_Unknowns(const DATA<TV>& data,const FORCE<TV>& force) const
{
    Matrix<T,Dynamic,1> velocities,forces,unknowns;
    data.Pack_Velocities(velocities);
    STORED_FORCE<T> stored_force;
    force.Pack_Forces(stored_force);
    forces=stored_force.Vector();
    unknowns.resize(velocities.size()+forces.size());
    for(int i=0;i<velocities.size();i++){
        unknowns(i)=velocities(i);}
    for(int i=0;i<forces.size();i++){
        unknowns(i+velocities.size())=forces(i);}
    return unknowns;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Increment_Unknowns(const Matrix<T,Dynamic,1>& unknowns,DATA<TV>& data,FORCE<TV>& force)
{
    int velocity_dof=Velocity_DOF();
    Matrix<T,Dynamic,1> solve_velocities=unknowns.block(0,0,velocity_dof,1);
    STORED_FORCE<T> stored_force;
    force.Pack_Forces(stored_force);
    stored_force.Set(unknowns.block(velocity_dof,0,unknowns.size()-velocity_dof,1));
    force.Increment_Forces(stored_force,1);

    Matrix<T,Dynamic,1> current_velocities;
    data.Pack_Velocities(current_velocities);
    Matrix<T,Dynamic,1> candidate_velocities=current_velocities+solve_velocities;
    Unpack_Velocities(data,candidate_velocities);
    data.Step();
}
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Identify_DOF(const DATA<TV>& data,const FORCE<TV>& force,int index)
{
    int current_index=0;
    for(int i=0;i<data.size();i++){
        int data_size=kinematic_projection_matrices[i].rows();//=data[i]->Velocity_DOF();
        if(index<current_index+data_size){
            data[i]->Identify_DOF(index-current_index);
            return;}
        current_index+=data_size;
    }
    for(int i=0;i<force.size();i++){
        int force_size=force[i]->DOF();
        if(index<current_index+force_size){
            force[i]->Identify_DOF(index-current_index);
            return;}
        current_index+=force_size;}
    LOG::cout<<"Unidentified DOF"<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Unpack_Velocities(DATA<TV>& data,const Matrix<T,Dynamic,1>& velocities)
{
    assert(data.size()==1);
    Matrix<Matrix<T,Dynamic,1>,Dynamic,1> full_velocities(data.size());
    int current_index=0;
    for(int i=0;i<data.size();i++){
        int compressed_size=kinematic_projection_matrices[i].rows();
        full_velocities[i]=kinematic_projection_matrices[i].transpose()*velocities.block(current_index,0,current_index+compressed_size,1);
        current_index+=compressed_size;
    }
    Matrix<T,Dynamic,1> full_velocity;
    Merge_Block_Vectors(full_velocities,full_velocity);
    data.Unpack_Velocities(full_velocity);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Store_Errors(DATA<TV>& data,const Matrix<T,Dynamic,1>& errors)
{
    assert(data.size()==1);
    Matrix<Matrix<T,Dynamic,1>,Dynamic,1> full_errors(data.size());
    int current_index=0;
    for(int i=0;i<data.size();i++){
        int compressed_size=kinematic_projection_matrices[i].rows();
        full_errors[i]=kinematic_projection_matrices[i].transpose()*errors.block(current_index,0,current_index+compressed_size,1);
        current_index+=compressed_size;
    }
    Matrix<T,Dynamic,1> full_error;
    Merge_Block_Vectors(full_errors,full_error);
    data.Store_Errors(full_error);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Initialize(DATA<TV>& data,FORCE<TV>& force)
{
    kinematic_projection_matrices.resize(data.size());
    for(int i=0;i<data.size();i++){    
        data[i]->Kinematic_Projection(kinematic_projection_matrices[i]);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void NONLINEAR_EQUATION<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic)
{
    system.Initialize(data,force);
    for(int i=0;i<data.size();i++){
        data[i]->Inertia(dt,system.jacobian_block_terms[i],system.inverse_inertia_matrices[i],system.error_blocks[i]);}
    for(int i=0;i<force.size();i++){
        force[i]->Identify_Interactions_And_Compute_Errors(data,force,dt,time,system,stochastic);}
    system.Scale_Errors(data,force,kinematic_projection_matrices);
    for(int i=0;i<force.size();i++){
        force[i]->Compute_Derivatives(data,force,system);}
    system.Scale_Derivatives(data,force,kinematic_projection_matrices);
    jacobian.resize(data.Velocity_DOF(),data.Velocity_DOF());
    Merge_Block_Matrices(system.jacobian_blocks,jacobian);
    Merge_Block_Vectors(system.error_blocks,error);

    // build the Hessian
    SparseMatrix<T> hessian_addition;
    Merge_Block_Matrices(system.hessian_blocks,hessian_addition);
    hessian=jacobian.adjoint()*jacobian+hessian_addition;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(NONLINEAR_EQUATION)
