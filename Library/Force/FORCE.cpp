#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> int FORCE<TV>::
Force_DOF() const
{
    int total_size=0;
    for(auto force_type : (*this)){total_size+=force_type->Force_DOF();}
    return total_size;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Pack_Forces(Matrix<T,Dynamic,1>& forces) const
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
Increment_Forces(const Matrix<T,Dynamic,1>& forces)
{
    int current_position=0;
    for(auto force_type : (*this)){
        int force_size=force_type->Force_DOF();
        force_type->Increment_Forces(forces.block(current_position,0,force_size,1));
        current_position+=force_size;
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void FORCE<TV>::
Viewer(const DATA<TV>& data,osg::Group*& root)
{
    for(auto force_type : (*this)){force_type->Viewer(data,root);}
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(FORCE)
