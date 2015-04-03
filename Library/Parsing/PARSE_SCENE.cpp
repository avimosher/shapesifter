#include <Parsing/PARSE_SCENE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/LOG.h>
#include <json/json.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> bool PARSE_SCENE<TV>::
Parse_Scene(std::istream& input,SIMULATION<TV>& simulation)
{
    Json::Value root;
    Json::Reader reader;
    bool success=reader.parse(input,root);
    LOG::cout<<"Parse success: "<<success<<std::endl;
    if(!success){
        LOG::cout<<reader.getFormattedErrorMessages()<<std::endl;
        return false;
    }
    Json::Value arrayRoot=root["root"];
    Parse_String(root["output_directory"],simulation.output_directory);
    Parse_Scalar(root["dt"],simulation.dt);
    Parse_Scalar(root["last_time"],simulation.last_time);
    for(int i=0;i<arrayRoot.size();i++){PARSER_REGISTRY<TV,void>::Parse(arrayRoot[i],simulation);}
    return true;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(PARSE_SCENE)
