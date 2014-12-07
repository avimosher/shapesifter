///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class DATA
/////////////////////////////////////////////////////////////////////// 
#include <Data/DATA.h>
#include <iostream>

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
template<class TV> void DATA<TV>::
Variables(Matrix<typename TV::Scalar,Dynamic,1>& variables)
{
    int total_size=0;
    for (DATA_TYPE<TV>* data_type : (*this)){
        total_size+=data_type->Size();
    }
    variables.resize(total_size,1);
    int current_position=0;
    for (DATA_TYPE<TV>* data_type : (*this)){
        int data_size=data_type->Size();
        variables.block(current_position,0,data_size,1)=data_type->Variables();
        current_position+=data_size;
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Step(QUALITY& step_quality,Matrix<T,Dynamic,1> solve_result)
{
    int current_position=0;
    for(DATA_TYPE<TV>* data_type : (*this)){
        int data_size=data_type->Size();
        std::cout<<"Stepcaller"<<std::endl;
        data_type->Step(solve_result.block(current_position,0,data_size,1));
        current_position+=data_size;
    }
}
/////////////////////////////////////////////////////////////////////// 
GENERIC_TYPE_DEFINITION(DATA)
