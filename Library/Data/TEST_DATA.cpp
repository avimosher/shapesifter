#include <Data/DATA.h>
#include <Data/TEST_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(TEST_DATA)
GENERIC_CEREAL_REGISTRATION(TEST_DATA)
DEFINE_AND_REGISTER_PARSER(TEST_DATA,void)
{
    auto test_data=std::make_shared<TEST_DATA<TV>>();
    test_data->internal_data=node["internal_data"].asDouble();
    std::cout<<test_data->internal_data<<std::endl;
    simulation.data.insert({test_data->Name(),test_data});
    return 0;
}
