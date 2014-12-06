//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class SIMPLE_REPULSION_FORCE
///////////////////////////////////////////////////////////////////////
#include <Force/SIMPLE_REPULSION_FORCE.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> SIMPLE_REPULSION_FORCE<TV>::
SIMPLE_REPULSION_FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> SIMPLE_REPULSION_FORCE<TV>::
~SIMPLE_REPULSION_FORCE()
{}
///////////////////////////////////////////////////////////////////////
template<class TV> void SIMPLE_REPULSION_FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side)
{
#if 0
    FRAME<TV> frame1=data.Updated_Frame(rigid_structure1->frame,rigid_velocity->twist(body_index1));
    FRAME<TV> frame2=data.Updated_Frame(rigid_structure2->frame,rigid_velocity->twist(body_index2));
    TV direction=frame1.t-frame2.t;
    T distance=direction.Normalize();
    for(int axis=0;axis<TV::RowsAtCompileTime;axis++) {
        projection(i)(1,axis)=direction(axis);
    }
    right_hand_side(i)=constraint.target_distance-distance;
    C_linear=projection*translation_map;
#endif
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SIMPLE_REPULSION_FORCE)
