#include <Driver/SIMULATION.h>
#include <Evolution/EVOLUTION.h>
#include <Evolution/EVOLUTION_STEP.h>
#include <Parsing/PARSE_NONLINEAR_EQUATION.h>
#include <Parsing/PARSER_REGISTRY.h>

using namespace Mechanics;
template<class TV> void PARSE_NONLINEAR_EQUATION<TV>::
Parse(Json::Value& node,SIMULATION<TV>& simulation)
{
    auto evolution_step=std::make_shared<EVOLUTION_STEP<TV>>();
    evolution_step->equation=new NONLINEAR_EQUATION<TV>();
    simulation.evolution.push_back(evolution_step);
}

REGISTER_PARSER(PARSE_NONLINEAR_EQUATION)
GENERIC_TYPE_DEFINITION(PARSE_NONLINEAR_EQUATION)
