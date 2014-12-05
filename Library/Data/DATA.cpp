///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class DATA
/////////////////////////////////////////////////////////////////////// 
#include <Data/DATA.h>
using namespace Mechanics;
template<class TV> DATA<TV>::
DATA()
{
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> DATA<TV>::
Variables()
{
    int total_size=0;
    std::vector<Matrix<T,Dynamic,1>> data_variables; //TODO: preallocate
    for (DATA_TYPE<TV>* data_type : data){
        data_variables.push_back(data_type->Variables());
        total_size+=data_variables.back().rows();
    }
    Matrix<T,Dynamic,1> variables(total_size);
    int current_position=0;
    for(Matrix<T,Dynamic,1> data_variable : data_variables){
        variables.block(current_position,current_position+data_variable.rows(),0,0)=data_variable;
        current_position+=data_variable.rows();
    }
    return variables;
}
/////////////////////////////////////////////////////////////////////// 
GENERIC_TYPE_DEFINITION(DATA)
