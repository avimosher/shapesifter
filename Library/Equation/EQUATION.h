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

    virtual void Unpack_Velocities(DATA<TV>& data,const Matrix<T,Dynamic,1>& velocities){};
    virtual void Store_Errors(DATA<TV>& data,const Matrix<T,Dynamic,1>& errors){};
    virtual void Initialize(DATA<TV>& data,FORCE<TV>& force){};
    virtual void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic)=0;
    virtual void Linearize_Around(const Vector& x){};
    virtual T Evaluate()=0;
    virtual void Gradient(Matrix<T,Dynamic,1>& gradient) const=0;
    virtual void RHS(Matrix<T,Dynamic,1>& rhs) const=0;
    virtual void Hessian(SparseMatrix<T>& hessian) const=0;
    virtual void Jacobian(SparseMatrix<T>& jacobian) const=0;
    virtual int System_Size()=0;
    virtual int Velocity_DOF() const=0;
};
}
#endif
