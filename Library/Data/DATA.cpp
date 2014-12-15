///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class DATA
/////////////////////////////////////////////////////////////////////// 
#include <Data/DATA.h>
#include <iostream>
#include <fstream>
#include <cereal/archives/json.hpp>
#include <osg/Group>

using namespace Mechanics;
/////////////////////////////////////////////////////////////////////// 
template<class TV> DATA<TV>::
DATA()
    :output_directory(".")
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
    for (std::unique_ptr<DATA_TYPE<TV>>& data_type : (*this)){
        total_size+=data_type->Size();
    }
    variables.resize(total_size,1);
    int current_position=0;
    for (std::unique_ptr<DATA_TYPE<TV>>& data_type : (*this)){
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
    for(std::unique_ptr<DATA_TYPE<TV>>& data_type : (*this)){
        int data_size=data_type->Size();
        std::cout<<"Stepcaller"<<std::endl;
        data_type->Step(solve_result.block(current_position,0,data_size,1));
        current_position+=data_size;
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Write(const int frame)
{
    std::ostringstream stringStream;
    stringStream<<output_directory<<"/frame."<<frame;
    std::ofstream output(stringStream.str().c_str(),std::ios::out);
    cereal::JSONOutputArchive archive(output);
    archive((*this)[0]);
    //archive(*((std::vector<std::unique_ptr<DATA_TYPE<TV>>>*)this));
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void DATA<TV>::
Read(const int frame)
{
    std::ostringstream stringStream;
    stringStream<<output_directory<<"/frame."<<frame;
    std::ifstream input(stringStream.str().c_str(),std::ios::in);
    cereal::JSONInputArchive archive(input);
    archive((*this)[0]);
    //archive(*((std::vector<std::unique_ptr<DATA_TYPE<TV>>>*)this));
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
    root=root?root:new osg::Group();
    for(std::unique_ptr<DATA_TYPE<TV>>& data_type : (*this)){
        data_type->Viewer(root);
    }
}
/////////////////////////////////////////////////////////////////////// 
GENERIC_TYPE_DEFINITION(DATA)
