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

    virtual T Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,const bool stochastic)=0;
    virtual Matrix<T,Dynamic,1> Solve()=0;
    virtual T Sufficient_Descent_Factor(const Matrix<T,Dynamic,1>& direction)=0;
};
}
#endif
