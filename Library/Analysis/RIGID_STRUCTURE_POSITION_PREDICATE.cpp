#include <Analysis/RIGID_STRUCTURE_POSITION_PREDICATE.h>
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> TV RIGID_STRUCTURE_POSITION_PREDICATE<TV>::
Vector(const SIMULATION<TV>& simulation)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto rigid_structure=rigid_data->Structure(name);
    return rigid_structure->frame*offset;
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_POSITION_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(RIGID_STRUCTURE_POSITION_PREDICATE,PREDICATE)
{
    auto predicate=std::make_shared<RIGID_STRUCTURE_POSITION_PREDICATE<TV>>();
    predicate->name=node["name"].asString();
    Parse_Vector(node["offset"],predicate->offset,TV());
    return predicate;
}
