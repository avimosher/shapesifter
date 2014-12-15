#include <Parsing/FACTORY.h>
#include <Parsing/PARSE_SCENE.h>
#include <json/json.h>

using namespace Mechanics;

template<class TV> bool PARSE_SCENE<TV>::
Parse_Scene(std::istream& input,DATA<TV>& data)
{
    Json::Value root;
    Json::Reader reader;
    bool success=reader.parse(input,root);
    std::cout<<"Parse success: "<<success<<std::endl;
    if(!success){return false;}
    Json::Value arrayRoot=root["root"];
    std::cout<<arrayRoot.size()<<std::endl;
    for(int i=0;i<arrayRoot.size();i++){FACTORY<TV>::Parse(arrayRoot[i],data);}
    return true;
}

GENERIC_TYPE_DEFINITION(PARSE_SCENE)
