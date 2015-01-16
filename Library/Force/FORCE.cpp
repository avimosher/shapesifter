//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class FORCE
///////////////////////////////////////////////////////////////////////
#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> int FORCE<TV>::
Force_DOF() const
{
    int total_size=0;
    for(auto force_type : (*this)){
        total_size+=force_type->Force_DOF();
    }
    return total_size;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Pack_Forces(Matrix<T,Dynamic,1>& forces)
{
    forces.resize(Force_DOF(),1);
    int current_position=0;
    for(auto force_type : (*this)){
        int force_size=force_type->Force_DOF();
        Block<Matrix<T,Dynamic,1>> block=forces.block(current_position,0,force_size,1);
        force_type->Pack_Forces(block);
        current_position+=force_size;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Unpack_Forces(const Matrix<T,Dynamic,1>& forces)
{
    int current_position=0;
    for(auto force_type : (*this)){
        int force_size=force_type->Force_DOF();
        force_type->Unpack_Forces(forces.block(current_position,0,force_size,1));
        current_position+=force_size;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side)
{
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(FORCE)
