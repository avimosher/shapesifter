//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class FORCE
//#####################################################################
#ifndef __FORCE__
#define __FORCE__

#include <Utilities/TYPE_UTILITIES.h>
#include <memory>
#include <vector>
#include <Eigen/SparseCore>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE_TYPE;


template<class TV>
class FORCE:public std::vector<std::shared_ptr<FORCE_TYPE<TV>>>
{
    typedef typename TV::Scalar T;

public:
    FORCE();
    ~FORCE();

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    int Force_DOF() const;
    void Pack_Forces(Matrix<T,Dynamic,1>& forces);
    void Unpack_Forces(const Matrix<T,Dynamic,1>& forces);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side);
};
}
#endif
