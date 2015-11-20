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
    inverse_inertia_matrices.resize(data.size());
    for(int i=0;i<data.size();i++){
        data[i]->Inertia(dt,system.matrix_block_terms[i],inverse_inertia_matrices[i],system.right_hand_side_blocks[i]);}
    for(int i=0;i<force.size();i++){
        force[i]->Linearize(data,force,dt,time,system,stochastic);}
    system.Scale_Blocks(data,force,kinematic_projection_matrices,inverse_inertia_matrices);
    jacobian.resize(data.Velocity_DOF(),data.Velocity_DOF());
    Merge_Block_Matrices(system.jacobian_blocks,jacobian);
    Merge_Block_Vectors(system.right_hand_side_blocks,right_hand_side);

    // build the Hessian
    //hessian=jacobian.adjoint()*jacobian;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(NONLINEAR_EQUATION)
