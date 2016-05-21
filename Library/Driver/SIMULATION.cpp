#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Evolution/EVOLUTION.h>
#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <fstream>
#include <Eigen/Geometry>
#include <sys/stat.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> SIMULATION<TV>::
SIMULATION()
    :data(*new DATA<TV>()),evolution(*new EVOLUTION<TV>()),force(*new FORCE<TV>()),
    current_frame(0),restart_frame(0),output_number(0),dt(.1),time(0),restart(false),substeps(false),write(true),
    output_directory("."),title("")
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
Write(const std::string& frame_title)
{
    if(write){
        title=frame_title;
        std::ostringstream stringStream;
        mkdir(output_directory.c_str(),0777);
        stringStream<<output_directory<<"/frame."<<output_number++;
        std::ofstream output(stringStream.str().c_str(),std::ios::out);
        //cereal::BinaryOutputArchive archive(output);
        cereal::JSONOutputArchive archive(output);
        archive(time,title,data,data.random,force);
    }
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> bool SIMULATION<TV>::
Read(const int frame)
{
    std::ostringstream stringStream;
    stringStream<<output_directory<<"/frame."<<frame;
    std::ifstream input(stringStream.str().c_str(),std::ios::in);
    if(!input.is_open()){return false;}
    //cereal::BinaryInputArchive archive(input);
    cereal::JSONInputArchive archive(input);
    archive(time,title,data,data.random,force);
    return true;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SIMULATION)
