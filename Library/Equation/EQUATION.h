#ifndef __EQUATION__
#define __EQUATION__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE;
template<class T> class QUALITY;

template<class TV>
class EQUATION
{
    typedef typename TV::Scalar T;
    typedef Matrix<T,Dynamic,1> Vector;
public:
    EQUATION(){};
    ~EQUATION(){};

    virtual void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic)=0;
    virtual void Linearize_Around(const Vector& x){};
    virtual T Evaluate()=0;
    virtual Matrix<T,Dynamic,1> Gradient()=0;
    virtual Matrix<T,Dynamic,1> RHS()=0;
    virtual SparseMatrix<T> Hessian()=0;
    virtual SparseMatrix<T> Jacobian()=0;
    virtual int System_Size()=0;
};
}
#endif
