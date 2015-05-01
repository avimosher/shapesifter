#include <Data/DATA.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/RANDOM.h>
#include <fstream>
#include <iostream>
#include <osg/Group>
using namespace Mechanics;
/////////////////////////////////////////////////////////////////////// 
template<class TV> DATA<TV>::
DATA()
    :random(*new RANDOM<T>)
{
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> DATA<TV>::
~DATA()
{
    delete &random;
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> int DATA<TV>::
Velocity_DOF() const
{
    int total_size=0;
    for(auto data_type : (*this)){total_size+=data_type->Velocity_DOF();}
    return total_size;
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> int DATA<TV>::
Position_DOF() const
{
    int total_size=0;
    for(auto data_type : (*this)){total_size+=data_type->Position_DOF();}
    return total_size;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DATA<TV>::
Pack_Velocities(Matrix<T,Dynamic,1>& velocities)
{
    velocities.resize(Velocity_DOF());
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type->Velocity_DOF();
        Block<Matrix<T,Dynamic,1>> block=velocities.block(current_position,0,data_size,1);
        data_type->Pack_Velocities(block);
        current_position+=data_size;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DATA<TV>::
Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities)
{
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type->Velocity_DOF();
        data_type->Unpack_Velocities(velocities.block(current_position,0,data_size,1));
        current_position+=data_size;}
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Pack_Positions(Matrix<T,Dynamic,1>& positions)
{
    positions.resize(Position_DOF());
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type->Position_DOF();
        Block<Matrix<T,Dynamic,1>> block=positions.block(current_position,0,data_size,1);
        data_type->Pack_Positions(block);
        current_position+=data_size;}
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Unpack_Positions(const Matrix<T,Dynamic,1>& positions)
{
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type->Position_DOF();
        data_type->Unpack_Positions(positions.block(current_position,0,data_size,1));
        current_position+=data_size;}
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Step()
{
    for(auto data_type : (*this)){data_type->Step(*this);}
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Viewer(osg::Group*& root)
{
    for(auto data_type : (*this)){data_type->Viewer(root);}
}
/////////////////////////////////////////////////////////////////////// 
GENERIC_TYPE_DEFINITION(DATA)
DEFINE_AND_REGISTER_PARSER(DATA,void)
{
    for(Json::ValueIterator itr=node.begin();itr!=node.end();itr++){
        if(itr.key()=="type"){continue;}
        if(itr.key()=="domain"){
            simulation.data.periodic=true;
            Parse_Vector(node["domain"]["minimum_corner"],simulation.data.domain.minimum_corner);
            Parse_Vector(node["domain"]["maximum_corner"],simulation.data.domain.maximum_corner);
            continue;}
        if(itr.key()!="type"){Parse_Scalar(*itr,simulation.data.globals[itr.key().asString()]);}}
    return 0;
}
