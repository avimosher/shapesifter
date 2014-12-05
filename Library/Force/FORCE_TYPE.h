//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class FORCE_TYPE
//#####################################################################
#ifndef __FORCE_TYPE__
#define __FORCE_TYPE__

#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/SparseCore>

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

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side);
};
}
#endif
