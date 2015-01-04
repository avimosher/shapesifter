//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class FORCE_TYPE
//#####################################################################
#ifndef __FORCE_TYPE__
#define __FORCE_TYPE__

#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/SparseCore>
#include <limits>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE;

template<class TV>
class FORCE_TYPE
{
    typedef typename TV::Scalar T;

public:
    FORCE_TYPE();
    ~FORCE_TYPE();

    virtual T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time){return std::numeric_limits<T>::max();}
    virtual void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)=0;
};
}
#endif
