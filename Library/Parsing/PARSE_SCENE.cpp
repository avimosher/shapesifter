#include <Parsing/PARSE_SCENE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <json/json.h>

using namespace Mechanics;

template<class TV> bool PARSE_SCENE<TV>::
Parse_Scene(std::istream& input,SIMULATION<TV>& simulation)
{
    Json::Value root;
    Json::Reader reader;
    bool success=reader.parse(input,root);
    std::cout<<"Parse success: "<<success<<std::endl;
    if(!success){
        std::cout<<reader.getFormattedErrorMessages()<<std::endl;
        return false;
    }
    Json::Value arrayRoot=root["root"];
    std::cout<<arrayRoot.size()<<std::endl;
    for(int i=0;i<arrayRoot.size();i++){PARSER_REGISTRY<TV>::Parse(arrayRoot[i],simulation);}
    return true;
}

GENERIC_TYPE_DEFINITION(PARSE_SCENE)
