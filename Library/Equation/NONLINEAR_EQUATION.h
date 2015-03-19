#ifndef __NONLINEAR_EQUATION__
#define __NONLINEAR_EQUATION__

#include <Equation/EQUATION.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

namespace Mechanics{
template<class TV> class DATA;

template<class T>
class RowPreconditioner:public DiagonalPreconditioner<T>
{
public:
    using DiagonalPreconditioner<T>::m_invdiag;

    void SetDiagonal(const Matrix<T,Dynamic,1>& diagonal)
    {m_invdiag=diagonal.cwiseInverse();}
};

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
    Matrix<T,Dynamic,1> conditioner;

    NONLINEAR_EQUATION(){};
    ~NONLINEAR_EQUATION(){};

    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic);
    Matrix<T,Dynamic,1> Solve();
    Matrix<T,Dynamic,1> Solve_Trust_Region();
    T Sufficient_Descent_Factor(const Matrix<T,Dynamic,1>& direction);
    T Evaluate();
    Matrix<T,Dynamic,1> Gradient();
    SparseMatrix<T> Hessian();
    int System_Size();
    DEFINE_TYPE_NAME("NONLINEAR_EQUATION")
};
}
#endif
