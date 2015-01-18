///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class DATA
/////////////////////////////////////////////////////////////////////// 
#include <Data/DATA.h>
#include <iostream>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <osg/Group>

using namespace Mechanics;
/////////////////////////////////////////////////////////////////////// 
template<class TV> DATA<TV>::
DATA()
{
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> DATA<TV>::
~DATA()
{
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> int DATA<TV>::
Velocity_DOF() const
{
    int total_size=0;
    for(auto data_type : (*this)){
        total_size+=data_type.second->Velocity_DOF();
    }
    return total_size;
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Variables(Matrix<typename TV::Scalar,Dynamic,1>& variables)
{
    int total_size=0;
    for(auto data_type : (*this)){
        total_size+=data_type.second->Size();
    }
    variables.resize(total_size,1);
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type.second->Size();
        variables.block(current_position,0,data_size,1)=data_type.second->Variables();
        current_position+=data_size;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DATA<TV>::
Pack_Velocities(Matrix<T,Dynamic,1>& velocities)
{
    velocities.resize(Velocity_DOF(),1);
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type.second->Velocity_DOF();
        Block<Matrix<T,Dynamic,1>> block=velocities.block(current_position,0,data_size,1);
        data_type.second->Pack_Velocities(block);
        current_position+=data_size;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void DATA<TV>::
Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities)
{
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type.second->Velocity_DOF();
        data_type.second->Unpack_Velocities(velocities.block(current_position,0,data_size,1));
        current_position+=data_size;
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Pack_Positions(Matrix<T,Dynamic,1>& positions)
{
    int total_size=0;
    for(auto data_type : (*this)){
        total_size+=data_type.second->Position_DOF();
    }
    positions.resize(total_size,1);
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type.second->Position_DOF();
        Block<Matrix<T,Dynamic,1>> block=positions.block(current_position,0,data_size,1);
        data_type.second->Pack_Positions(block);
        current_position+=data_size;
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Unpack_Positions(const Matrix<T,Dynamic,1>& positions)
{
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type.second->Position_DOF();
        data_type.second->Unpack_Positions(positions.block(current_position,0,data_size,1));
        current_position+=data_size;
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Step(QUALITY& step_quality,Matrix<T,Dynamic,1> solve_result)
{
    int current_position=0;
    for(auto data_type : (*this)){
        int data_size=data_type.second->Size();
        data_type.second->Step(*this,solve_result.block(current_position,0,data_size,1));
        current_position+=data_size;
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> typename TV::Scalar DATA<TV>::
Print_All()
{
    return (*this)[0]->Print();
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Viewer(osg::Group*& root)
{
    for(auto data_type : (*this)){
        data_type.second->Viewer(root);
    }
}
/////////////////////////////////////////////////////////////////////// 
GENERIC_TYPE_DEFINITION(DATA)
