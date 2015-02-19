#ifndef __NONLINEAR_EQUATION__
#define __NONLINEAR_EQUATION__

#include <Equation/EQUATION.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/Sparse>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class NONLINEAR_EQUATION : public EQUATION<TV>
{
    typedef typename TV::Scalar T;
public:
    SparseMatrix<T> matrix;
    Matrix<SparseMatrix<T>,Dynamic,Dynamic> full_matrix;
    Matrix<T,Dynamic,1> right_hand_side;
    Matrix<T,Dynamic,1> right_hand_side_full;
    SparseMatrix<T> J;
    Matrix<Matrix<T,Dynamic,1>,Dynamic,1> full_right_hand_side;

    NONLINEAR_EQUATION(){};
    ~NONLINEAR_EQUATION(){};

    T Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic);
    Matrix<T,Dynamic,1> Solve();
    T Sufficient_Descent_Factor(const Matrix<T,Dynamic,1>& direction);
    DEFINE_TYPE_NAME("NONLINEAR_EQUATION")
};
}
#endif
