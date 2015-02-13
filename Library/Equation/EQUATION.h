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

public:
    EQUATION(){};
    ~EQUATION(){};

    virtual void Linearize(DATA<TV>& data,FORCE<TV>& force,const Matrix<T,Dynamic,1>& velocities,const T dt,const T time,const bool stochastic)=0;
    virtual T Calculate_RHS_And_Norm(const DATA<TV>& data,const FORCE<TV>& force,const Matrix<T,Dynamic,1>& velocities)=0;
    virtual Matrix<T,Dynamic,1> Solve(const Matrix<T,Dynamic,1>& guess)=0;
    virtual bool Satisfied(DATA<TV>& data,FORCE<TV>& force,const Matrix<T,Dynamic,1>& solve_result,QUALITY<T>& solve_quality)=0;
};
}
#endif
