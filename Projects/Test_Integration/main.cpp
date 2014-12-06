#include <Data/DATA.h>
#include <Driver/DRIVER.h>
#include <Evolution/EVOLUTION.h>
#include <Force/FORCE.h>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <string>

using namespace Mechanics;
int main()
{
    typedef double T;
    typedef Matrix<T,1,1> TV;
    DATA<TV> data;
    EVOLUTION<TV> evolution;
    FORCE<TV> force;
    DRIVER<TV> driver(data,evolution,force);
    std::cout<<driver.restart_frame<<std::endl;
    return 0;
}
