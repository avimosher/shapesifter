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
    std::vector<Matrix<T,Dynamic,1>> data_variables; //TODO: preallocate
    for (DATA_TYPE<TV>* data_type : (*this)){
        data_variables.push_back(data_type->Variables());
        total_size+=data_variables.back().rows();
    }
    variables.resize(total_size,1);
    std::cout<<total_size<<std::endl;
    int current_position=0;
    for(Matrix<T,Dynamic,1> data_variable : data_variables){
        variables.block(current_position,0,data_variable.rows(),1)=data_variable;
        current_position+=data_variable.rows();
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Step(QUALITY& step_quality,Matrix<T,Dynamic,1> solve_result)
{
}
/////////////////////////////////////////////////////////////////////// 
GENERIC_TYPE_DEFINITION(DATA)
