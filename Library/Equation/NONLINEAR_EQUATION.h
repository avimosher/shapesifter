#ifndef __NONLINEAR_EQUATION__
#define __NONLINEAR_EQUATION__

#include <Equation/EQUATION.h>
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
    SparseMatrix<T> matrix;
    Matrix<SparseMatrix<T>,Dynamic,Dynamic> inverse_inertia_matrix;
    Matrix<SparseMatrix<T>,Dynamic,Dynamic> full_matrix;
    Matrix<T,Dynamic,1> right_hand_side;
    SparseMatrix<T> jacobian;
    SparseMatrix<T> inverse_inertia;
    Matrix<Matrix<T,Dynamic,1>,Dynamic,1> full_right_hand_side;

    NONLINEAR_EQUATION(){};
    ~NONLINEAR_EQUATION(){};

    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic);
    T Evaluate();
    Matrix<T,Dynamic,1> Gradient();
    SparseMatrix<T> Hessian();
    int System_Size();
    DEFINE_TYPE_NAME("NONLINEAR_EQUATION")
};
}
#endif
