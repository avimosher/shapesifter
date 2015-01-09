//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
template<class TV> RIGID_STRUCTURE<TV>::
RIGID_STRUCTURE()
{
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RIGID_STRUCTURE)
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE,void)
{
    std::shared_ptr<RIGID_STRUCTURE<TV>> structure=std::make_shared<RIGID_STRUCTURE<TV>>();
    Json::Value position=node["position"];
    for(int i=0;i<position.size();i++){structure->frame.position[i]=position[i].asDouble();}
    structure->name=node["name"].asString();
    auto data_element=simulation.data.find(RIGID_STRUCTURE_DATA<TV>::Static_Name());
    auto rigid_structure_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data_element->second);
    rigid_structure_data->structures.push_back(structure);
    return 0;
}
