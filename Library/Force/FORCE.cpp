//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class FORCE
///////////////////////////////////////////////////////////////////////
#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> FORCE<TV>::
FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> FORCE<TV>::
~FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side)
{
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(FORCE)