#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
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
        int data_size=data[i]->Velocity_DOF();
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
    int full_size=data.size()+force.size();
    full_matrix.resize(full_size,full_size);
    full_right_hand_side.resize(full_size,1);

    std::vector<std::vector<Triplet<T>>> force_terms;

    full_right_hand_side[0].resize(data.Velocity_DOF(),1);
    full_right_hand_side[0].setZero();
    int size=full_right_hand_side[0].rows();
    jacobian.resize(size,size);
    full_matrix(0,0).resize(size,size);
    // TODO: one block matrix per data type too
    force_terms.resize(data.size());
    inverse_inertia_matrices.resize(data.size());
    for(int i=0;i<data.size();i++){
        data[i]->Inertia(dt,force_terms[i],inverse_inertia_matrices[i],full_right_hand_side[i]);}

    for(int i=0;i<force.size();i++){
        // TODO: force_terms needs to be properly handled when there are multiple data types
        force[i]->Linearize(data,dt,time,force_terms[0],full_matrix(i+1,0),full_matrix(0,i+1),full_right_hand_side[0],full_right_hand_side[i+1],stochastic);}
    for(int i=0;i<data.size();i++){
        full_matrix(i,i).setFromTriplets(force_terms[i].begin(),force_terms[i].end());
        full_right_hand_side[i]=kinematic_projection_matrices[i]*inverse_inertia_matrices[i]*full_right_hand_side[i];
        for(int j=0;j<data.size()+force.size();j++){
            full_matrix(i,j)=inverse_inertia_matrices[i]*full_matrix(i,j);}}
    for(int i=0;i<data.size()+force.size();i++){
        for(int j=0;j<data.size()+force.size();j++){
            if(i<data.size() && j<data.size()){
                full_matrix(i,j)=kinematic_projection_matrices[i]*full_matrix(i,j)*kinematic_projection_matrices[j].transpose();}
            else if(i<data.size()){
                full_matrix(i,j)=kinematic_projection_matrices[i]*full_matrix(i,j);}
            else if(j<data.size()){
                full_matrix(i,j)=full_matrix(i,j)*kinematic_projection_matrices[j].transpose();}}}

    // scale the jacobian and rhs according to scaling factor on f.  Can't vary with x.  Make it 1/inertia diagonal for force balance, 1 for constraints
    Merge_Block_Matrices(full_matrix,jacobian);
    // TODO: this MAY BE able to scale just the mass matrix terms; the forces are free to be rescaled.  HOWEVER this may mess up the accumulation.
    Merge_Block_Vectors(full_right_hand_side,right_hand_side);
    //LOG::cout<<"Jacobian before inertia: "<<std::endl<<jacobian<<std::endl;
    //jacobian=inverse_inertia*jacobian;
    //LOG::cout<<"RHS before inertia: "<<std::endl<<right_hand_side.transpose()<<std::endl;
    //right_hand_side=inverse_inertia*right_hand_side;
    //LOG::cout<<"RHS: "<<std::endl<<right_hand_side.transpose()<<std::endl;
    /*int index;
    LOG::cout<<"Min value at index "<<index<<" is "<<right_hand_side.minCoeff(&index)<<std::endl;
    Identify_DOF(data,force,index);
    LOG::cout<<"Max value at index "<<index<<" is "<<right_hand_side.maxCoeff(&index)<<std::endl;
    Identify_DOF(data,force,index);*/
    //LOG::cout<<"Gradient: "<<(-jacobian.adjoint()*right_hand_side).transpose()<<std::endl;
    //LOG::cout<<"Jacobian: "<<std::endl<<jacobian<<std::endl;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(NONLINEAR_EQUATION)
