#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
#include <Force/STORED_FORCE.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Pack_Forces(STORED_FORCE<T>& stored_force) const
{
    if(stored_force.Size()!=(*this).size()){
        stored_force.Resize((*this).size());
        for(int i=0;i<(*this).size();i++){stored_force[i]=(*this)[i]->Create_Stored_Force();}}
    for(int i=0;i<this->size();i++){
        auto force_type=(*this)[i];
        force_type->Pack_Forces(stored_force[i]);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Unpack_Forces(const STORED_FORCE<T>& stored_force)
{
    for(int i=0;i<this->size();i++){
        auto force_type=(*this)[i];
        force_type->Unpack_Forces(stored_force[i]);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Increment_Forces(const STORED_FORCE<T>& stored_force,T ratio)
{
    for(int i=0;i<this->size();i++){
        auto force_type=(*this)[i];
        force_type->Increment_Forces(stored_force[i],ratio);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void FORCE<TV>::
Store_Errors(const STORED_FORCE<T>& stored_error)
{
    for(int i=0;i<this->size();i++){
        auto force_type=(*this)[i];
        force_type->Store_Errors(stored_error[i]);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> bool FORCE<TV>::
Equations_Changed() const
{
    for(auto force : (*this)){
        if(force->Equations_Changed()){return true;}}
    return false;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(FORCE)
