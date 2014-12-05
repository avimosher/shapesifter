//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class TEST_FORCE
///////////////////////////////////////////////////////////////////////
#include <Force/TEST_FORCE.h>
using namespace Mechanics;
template<class TV> TEST_FORCE<TV>::
TEST_FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TEST_FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side)
{
    // let's say the test force is F=-kv, and v is our only system variable (?)
    // so v=v*-k*dt*v
    T k=4;
    force_terms.push_back(Triplet<T>(0,0,-k*dt));
    //right_hand_side(0)=-
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(TEST_FORCE)
