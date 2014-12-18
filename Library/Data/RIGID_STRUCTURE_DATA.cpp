//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE_DATA
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> RIGID_STRUCTURE_DATA<TV>::
RIGID_STRUCTURE_DATA()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Size()
{
    return this->size()*RIGID_STRUCTURE<TV>::STATIC_SIZE;
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> RIGID_STRUCTURE_DATA<TV>::
Variables()
{
    static int blockSize=6;
    Matrix<T,Dynamic,1> packed(RIGID_STRUCTURE<TV>::STATIC_SIZE*this->size(),1);
    /*for(auto structure : (*this)) {
        packed.block<blockSize,1>()=structure.Pack();
        }*/
    return packed;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Step(const Matrix<T,Dynamic,1>& variables)
{
}
///////////////////////////////////////////////////////////////////////
#if 0
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> RIGID_STRUCTURE_DATA<TV>::
Pack()
{
    Matrix<T,Dynamic,1> packed(blockSize*structures.size(),1);
    for(RIGID_STRUCTURE<TV>* structure : structures) {
        packed.block<blockSize,1>()=structure.Pack();
    }
    return packed;
}
#endif
///////////////////////////////////////////////////////////////////////
#if 0
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> RIGID_STRUCTURE_DATA<TV>::
Unpack(Matrix<T,Dynamic,1>& packed)
{
    for(RIGID_STRUCTURE<TV>* structure : structures) {
        structure.Unpack(packed.block<blockSize,1>());
    }
}
#endif
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_DATA)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE_DATA)
{
    simulation.data.insert({RIGID_STRUCTURE_DATA<TV>::Static_Name(),std::make_shared<RIGID_STRUCTURE_DATA<TV>>()});
}
