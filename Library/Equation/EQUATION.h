//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class EQUATION
//#####################################################################
#ifndef __EQUATION__
#define __EQUATION__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE;

template<class TV>
class EQUATION
{
    typedef typename TV::Scalar T;

public:
    EQUATION(){};
    ~EQUATION(){};

    virtual void Linearize(DATA<TV>& data,FORCE<TV>& force,const Matrix<T,Dynamic,1>& velocities,const T dt,const T time)=0;
    virtual Matrix<T,Dynamic,1> Solve(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)=0;
    virtual bool Satisfied(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)=0;
};
}
#endif
