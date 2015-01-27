#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RIGID_STRUCTURE)
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE,void)
{
    std::shared_ptr<RIGID_STRUCTURE<TV>> structure=std::make_shared<RIGID_STRUCTURE<TV>>();
    Parse_Vector(node["position"],structure->frame.position,TV());
    Parse_Scalar(node["collision_extent"],structure->collision_extent);
    Parse_Scalar(node["radius"],structure->radius,(T)0);
    Parse_Scalar(node["collision_radius"],structure->collision_radius,structure->radius);
    structure->name=node["name"].asString();
    auto data_element=simulation.data.find(RIGID_STRUCTURE_DATA<TV>::Static_Name());
    auto rigid_structure_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data_element->second);
    rigid_structure_data->structures.push_back(structure);
    return 0;
}
