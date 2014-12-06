//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class FORCE
//#####################################################################
#ifndef __FORCE__
#define __FORCE__

#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/SparseCore>
#include <vector>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE_TYPE;


template<class TV>
class FORCE:public std::vector<FORCE_TYPE<TV>*>
{
    typedef typename TV::Scalar T;

public:
    FORCE();
    ~FORCE();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side);
};
}
#endif
