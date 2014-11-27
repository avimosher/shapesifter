//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class SIMPLE_REPULSION_FORCE
///////////////////////////////////////////////////////////////////////
#include <Data/SIMPLE_REPULSION_FORCE.h>
using namespace Mechanics;
template<class TV> SIMPLE_REPULSION_FORCE<TV>::
SIMPLE_REPULSION_FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> SIMPLE_REPULSION_FORCE<TV>::
Linearize()
{
    FRAME<TV> frame1=data.Updated_Frame(rigid_structure1->frame,rigid_velocity->twist(body_index1));
    FRAME<TV> frame2=data.Updated_Frame(rigid_structure2->frame,rigid_velocity->twist(body_index2));
    TV direction=frame1.t-frame2.t;
    T distance=direction.Normalize();
    for(int axis=0;axis<TV::RowsAtCompileTime;axis++) {
        projection(i)(1,axis)=direction(axis);
    }
    right_hand_side(i)=constraint.target_distance-distance;
    C_linear=projection*translation_map;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SIMPLE_REPULSION_FORCE)
