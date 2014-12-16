#include <Data/DATA.h>
#include <Data/TEST_DATA.h>
#include <Parsing/PARSE_TEST_DATA.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <iostream>

using namespace Mechanics;
REGISTER_PARSER(PARSE_TEST_DATA)

template<class TV> void PARSE_TEST_DATA<TV>::
Parse(Json::Value& node,DATA<TV>& data)
{
    auto test_data=std::make_shared<TEST_DATA<TV>>();
    test_data->internal_data=node["internal_data"].asDouble();
    std::cout<<test_data->internal_data<<std::endl;
    data.insert({test_data->Name(),test_data});
}

GENERIC_TYPE_DEFINITION(PARSE_TEST_DATA)
