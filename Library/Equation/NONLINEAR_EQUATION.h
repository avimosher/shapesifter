//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class NONLINEAR_EQUATION
//#####################################################################
#ifndef __NONLINEAR_EQUATION__
#define __NONLINEAR_EQUATION__

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class NONLINEAR_EQUATION
{
    typedef typename TV::Scalar T;

public:
    NONLINEAR_EQUATION();
    ~NONLINEAR_EQUATION();

    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    Matrix<T,Dynamic,1> Solve(DATA<TV>& data,FORCE<TV>& force,const T target_time);
};
}
#endif
