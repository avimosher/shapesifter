//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class BROWNIAN_FORCE
//#####################################################################
#ifndef __BROWNIAN_FORCE__
#define __BROWNIAN_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class BROWNIAN_FORCE : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;

public:
    BROWNIAN_FORCE();
    ~BROWNIAN_FORCE();

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side);    
};
}
#endif
