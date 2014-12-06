#include <Data/DATA.h>
#include <Data/TEST_DATA.h>
#include <Driver/DRIVER.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Evolution/EVOLUTION.h>
#include <Evolution/EVOLUTION_STEP.h>
#include <Force/FORCE.h>
#include <Force/TEST_FORCE.h>
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
    TEST_DATA<TV> test_data;
    EVOLUTION<TV> evolution;
    FORCE<TV> force;
    TEST_FORCE<TV> test_force;
    NONLINEAR_EQUATION<TV> equation;
    EVOLUTION_STEP<TV> step;
    step.equation=&equation;
    evolution.push_back(&step);
    force.push_back(&test_force);
    data.push_back(&test_data);
    DRIVER<TV> driver(data,evolution,force);
    driver.last_time=10;
    driver.Execute_Main_Program();
    return 0;
}
