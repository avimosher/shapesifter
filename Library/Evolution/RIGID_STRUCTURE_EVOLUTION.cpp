//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE_EVOLUTION
///////////////////////////////////////////////////////////////////////
#include <Data/RIGID_STRUCTURE_EVOLUTION.h>
using namespace Mechanics;
template<class TV> RIGID_STRUCTURE_EVOLUTION<TV>::
RIGID_STRUCTURE_EVOLUTION()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar RIGID_STRUCTURE_EVOLUTION<TV>::
Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T time,const T target_time)
{

}
///////////////////////////////////////////////////////////////////////
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_STRUCTURE<VECTOR<T,3> >& rigid_structure)
{
    typedef VECTOR<T,3> TV;ROTATION<TV>& R=rigid_structure.Rotation();const TV &w=rigid_structure.twist.angular,&L=rigid_structure.Angular_Momentum();
    TV rotate_amount=w-(T).5*dt*rigid_structure.World_Space_Inertia_Tensor_Inverse_Times(TV::Cross_Product(w,L));
    R=ROTATION<TV>::From_Rotation_Vector(dt*rotate_amount)*R;R.Normalize();
    rigid_structure.Update_Angular_Velocity(); // Note that the value of w changes here.
}
///////////////////////////////////////////////////////////////////////
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_STRUCTURE<VECTOR<T,2> >& rigid_structure)
{
    rigid_structure.Rotation()=ROTATION<VECTOR<T,2> >::From_Rotation_Vector(dt*rigid_structure.twist.angular)*rigid_structure.Rotation();rigid_structure.Rotation().Normalize();
}
///////////////////////////////////////////////////////////////////////
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_STRUCTURE<VECTOR<T,1> >& rigid_structure)
{}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_EVOLUTION<TV>::
Position_Step(DATA<TV>& data,FORCE<TV>& force,const T time,const T target_time)
{
    RIGID_STRUCTURE_DATA<TV>& rigid_structure_data=*data.template Find_Data<RIGID_STRUCTURE_DATA<TV>*>();

    for(RIGID_STRUCTURE<TV>& rigid_structure : rigid_structure_data.Structures()) {
        rigid_structure.frame.t=data.Wrap(rigid_structure.frame.t+dt*rigid_structure.twist.linear);
        Update_Rotation_Helper(dt,rigid_structure);        
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_EVOLUTION<TV>::
Velocity_Step(DATA<TV>& data,FORCE<TV>& force,const T time,const T target_time)
{
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_EVOLUTION)
