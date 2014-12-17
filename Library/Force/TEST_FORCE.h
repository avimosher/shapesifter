//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class TEST_FORCE
//#####################################################################
#ifndef __TEST_FORCE__
#define __TEST_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class TEST_FORCE : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;

public:
    TEST_FORCE();
    ~TEST_FORCE();

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side);
    DEFINE_TYPE_NAME("TEST_FORCE")
};
}
#endif
