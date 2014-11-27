//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE_DATA
///////////////////////////////////////////////////////////////////////
#include <Data/RIGID_STRUCTURE_DATA.h>
using namespace Mechanics;
template<class TV> RIGID_STRUCTURE_DATA<TV>::
RIGID_STRUCTURE_DATA()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> RIGID_STRUCTURE_DATA<TV>::
Pack()
{
    Matrix<T,Dynamic,1> packed(blockSize*structures.size(),1);
    for(RIGID_STRUCTURE<TV>* structure : structures) {
        packed.block<blockSize,1>()=structure.Pack();
    }
    return packed;
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> RIGID_STRUCTURE_DATA<TV>::
Unpack(Matrix<T,Dynamic,1>& packed)
{
    for(RIGID_STRUCTURE<TV>* structure : structures) {
        structure.Unpack(packed.block<blockSize,1>());
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_DATA)
