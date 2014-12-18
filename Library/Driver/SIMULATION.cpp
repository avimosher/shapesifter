///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class SIMULATION
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Evolution/EVOLUTION.h>
#include <Force/FORCE.h>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/string.hpp>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> SIMULATION<TV>::
SIMULATION()
    :data(*new DATA<TV>()),evolution(*new EVOLUTION<TV>()),force(*new FORCE<TV>()),
    current_frame(0),restart_frame(0),time(0),restart(false),output_directory(".")
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> SIMULATION<TV>::
~SIMULATION()
{
    delete &data;
    delete &evolution;
    delete &force;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SIMULATION<TV>::
Write(const int frame)
{
    std::ostringstream stringStream;
    stringStream<<output_directory<<"/frame."<<frame;
    std::ofstream output(stringStream.str().c_str(),std::ios::out);
    cereal::BinaryOutputArchive archive(output);
    archive(data);
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void SIMULATION<TV>::
Read(const int frame)
{
    std::ostringstream stringStream;
    stringStream<<output_directory<<"/frame."<<frame;
    std::ifstream input(stringStream.str().c_str(),std::ios::in);
    cereal::BinaryInputArchive archive(input);
    archive(data);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SIMULATION)
