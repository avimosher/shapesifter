#include <Driver/DRIVER.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSE_SCENE.h>
#include <algorithm>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <string>
using namespace Mechanics;

char* Get_Command_Option(char** begin,char** end,const std::string& option)
{
    char **itr=std::find(begin,end,option);
    if(itr!=end){
        if(++itr!=end){return *itr;}
        else{return *(--itr);}}
    return 0;
}

int main(int argc,char **argv)
{
    typedef double T;
    typedef Matrix<T,3,1> TV;
    auto simulation=std::make_shared<SIMULATION<TV>>();

    char *restart=Get_Command_Option(argv,argv+argc,"-restart");
    if(restart){simulation->Set_Restart(std::stoi(restart));}
    if(Get_Command_Option(argv,argv+argc,"-substeps")){
        simulation->substeps=true;
    }
    std::ifstream config(Get_Command_Option(argv,argv+argc,"-scene"),std::ifstream::in);
    if(!PARSE_SCENE<TV>::Parse_Scene(config,*simulation)){return 1;}

    DRIVER<TV> driver(simulation);
    simulation->last_time=10;
    driver.Execute_Main_Program();
    return 0;
}
