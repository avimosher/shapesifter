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
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic)
{
    int full_size=data.size()+force.size();
    inverse_inertia_matrix.resize(full_size,full_size);
    full_matrix.resize(full_size,full_size);
    full_right_hand_side.resize(full_size,1);

    std::vector<std::vector<Triplet<T>>> force_terms;
    std::vector<Triplet<T>> inverse_inertia_terms;
    std::vector<Triplet<T>> hessian_terms;

    full_right_hand_side[0].resize(data.Velocity_DOF(),1);
    full_right_hand_side[0].setZero();
    int size=full_right_hand_side[0].rows();
    jacobian.resize(size,size);
    full_matrix(0,0).resize(size,size);
    // TODO: one block matrix per data type too
    force_terms.resize(data.size());
    for(int i=0;i<data.size();i++){
        data[i]->Inertia(dt,force_terms[i],full_right_hand_side[i]);
        int data_size=data[i]->Velocity_DOF();
        for(int j=0;j<force_terms[i].size();j++){
            inverse_inertia_terms.push_back(Triplet<T>(force_terms[i][j].row(),force_terms[i][j].col(),1/fabs(force_terms[i][j].value())));
        }
    }
    int running_index=data[0]->Velocity_DOF();
    for(int i=0;i<force.size();i++){
        // TODO: force_terms needs to be properly handled when there are multiple data types
        force[i]->Linearize(data,dt,time,hessian_terms,force_terms[0],full_matrix(i+1,0),full_right_hand_side[0],full_right_hand_side[i+1],stochastic);
        full_matrix(0,i+1)=full_matrix(i+1,0).transpose();
        int force_size=full_right_hand_side[i+data.size()].size();
        for(int j=0;j<force_size;j++){
            inverse_inertia_terms.push_back(Triplet<T>(running_index+j,running_index+j,1));
        }
        running_index+=force_size;
    }
    for(int i=0;i<data.size();i++){
        full_matrix(i,i).setFromTriplets(force_terms[i].begin(),force_terms[i].end());
    }

    // scale the jacobian and rhs according to scaling factor on f.  Can't vary with x.  Make it 1/inertia diagonal for force balance, 1 for constraints
    
    Merge_Block_Matrices(full_matrix,jacobian);
    inverse_inertia.resize(running_index,running_index);
    inverse_inertia.setFromTriplets(inverse_inertia_terms.begin(),inverse_inertia_terms.end());
    Merge_Block_Vectors(full_right_hand_side,right_hand_side);
    jacobian=inverse_inertia*jacobian;
    //LOG::cout<<"RHS before inertia: "<<right_hand_side.transpose()<<std::endl;
    right_hand_side=inverse_inertia*right_hand_side;
    //LOG::cout<<"RHS: "<<right_hand_side.transpose()<<std::endl;
    //LOG::cout<<"Gradient: "<<(-jacobian.adjoint()*right_hand_side).transpose()<<std::endl;
    //LOG::cout<<"Jacobian: "<<std::endl<<jacobian<<std::endl;
    //LOG::cout<<"Jacobian adjoint: "<<std::endl<<jacobian.adjoint()<<std::endl;
    //LOG::cout<<"Rows: "<<jacobian.rows()<<" cols: "<<jacobian.cols()<<std::endl;
    //LOG::cout<<"RHS rows: "<<right_hand_side.rows()<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar NONLINEAR_EQUATION<TV>::
Evaluate()
{
    return right_hand_side.squaredNorm()/2;
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> NONLINEAR_EQUATION<TV>::
Gradient()
{
    return -jacobian.adjoint()*right_hand_side;
}
///////////////////////////////////////////////////////////////////////
template<class TV> SparseMatrix<typename TV::Scalar> NONLINEAR_EQUATION<TV>::
Hessian()
{
    return jacobian.adjoint()*jacobian;
}
///////////////////////////////////////////////////////////////////////
template<class TV> int NONLINEAR_EQUATION<TV>::
System_Size()
{
    return right_hand_side.size();
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(NONLINEAR_EQUATION)
