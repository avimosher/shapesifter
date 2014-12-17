#include <Driver/DRIVER.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSE_SCENE.h>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <string>

using namespace Mechanics;
int main()
{
    typedef double T;
    typedef Matrix<T,1,1> TV;
    auto simulation=std::make_shared<SIMULATION<TV>>();

    std::ifstream test_config("config.json",std::ifstream::in);
    PARSE_SCENE<TV>::Parse_Scene(test_config,*simulation);

    DRIVER<TV> driver(simulation);
    simulation->last_time=10;
    driver.Execute_Main_Program();
    return 0;
}
