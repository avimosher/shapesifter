#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/TEST_FORCE.h>
#include <Parsing/PARSE_TEST_FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>

using namespace Mechanics;
template<class TV> void PARSE_TEST_FORCE<TV>::
Parse(Json::Value& node,SIMULATION<TV>& simulation)
{
    auto test_force=std::make_shared<TEST_FORCE<TV>>();
    simulation.force.push_back(test_force);
}

REGISTER_PARSER(PARSE_TEST_FORCE)
GENERIC_TYPE_DEFINITION(PARSE_TEST_FORCE)
