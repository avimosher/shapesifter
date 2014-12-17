#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSE_TEST_DATA.h>
#include <Parsing/PARSER_REGISTRY.h>

using namespace Mechanics;
template<class TV> void PARSE_TEST_DATA<TV>::
Parse(Json::Value& node,SIMULATION<TV>& simulation)
{
    auto test_data=std::make_shared<TEST_DATA<TV>>();
    test_data->internal_data=node["internal_data"].asDouble();
    std::cout<<test_data->internal_data<<std::endl;
    simulation.data.insert({test_data->Name(),test_data});
}

REGISTER_PARSER(PARSE_TEST_DATA)
GENERIC_TYPE_DEFINITION(PARSE_TEST_DATA)
