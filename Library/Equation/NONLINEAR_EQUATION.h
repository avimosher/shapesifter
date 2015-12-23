#ifndef __NONLINEAR_EQUATION__
#define __NONLINEAR_EQUATION__

#include <Equation/EQUATION.h>
#include <Equation/MATRIX_BUNDLE.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class NONLINEAR_EQUATION : public EQUATION<TV>
{
    typedef typename TV::Scalar T;
public:
    MATRIX_BUNDLE<TV> system;

    std::vector<SparseMatrix<T>> inverse_inertia_matrices;
    std::vector<SparseMatrix<T>> kinematic_projection_matrices;
    Matrix<T,Dynamic,1> error;
    SparseMatrix<T> hessian;
    SparseMatrix<T> jacobian;
    SparseMatrix<T> inverse_inertia;

    NONLINEAR_EQUATION(){};
    ~NONLINEAR_EQUATION(){};

    int Velocity_DOF() const{
        int dof=0;
        for(int i=0;i<kinematic_projection_matrices.size();i++){
            dof+=kinematic_projection_matrices[i].rows();}
        return dof;
    }
    T Evaluate(){return error.squaredNorm()/2;}
    void Gradient(Matrix<T,Dynamic,1>& gradient) const{gradient=-jacobian.adjoint()*error;}
    void RHS(Matrix<T,Dynamic,1>& rhs) const{rhs=error;}
    void Hessian(SparseMatrix<T>& hessian_out) const{hessian_out=hessian;}
    void Jacobian(SparseMatrix<T>& jacobian_out) const{jacobian_out=jacobian;}
    int System_Size(){return error.size();}

    Matrix<T,Dynamic,1> Get_Unknowns(const DATA<TV>& data,const FORCE<TV>& force) const;
    void Increment_Unknowns(const Matrix<T,Dynamic,1>& unknowns,DATA<TV>& data,FORCE<TV>& force);
    void Identify_DOF(const DATA<TV>& data,const FORCE<TV>& force,int index);
    void Unpack_Velocities(DATA<TV>& data,const Matrix<T,Dynamic,1>& velocities);
    void Store_Errors(DATA<TV>& data,const Matrix<T,Dynamic,1>& errors);
    void Initialize(DATA<TV>& data,FORCE<TV>& force);
    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic);
    DEFINE_TYPE_NAME("NONLINEAR_EQUATION")
};
}
#endif
