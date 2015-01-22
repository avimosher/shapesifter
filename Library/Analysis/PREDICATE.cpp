#include <Analysis/PREDICATE.h>
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar DISTANCE_PREDICATE<TV>::
Scalar(const DATA<TV>& data,const FORCE<TV>& force)
{
    return data.Minimum_Offset(first_endpoint->Vector(data,force),second_endpoint->Vector(data,force)).norm();
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(DISTANCE_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(DISTANCE_PREDICATE,PREDICATE)
{
    auto predicate=std::make_shared<DISTANCE_PREDICATE<TV>>();
    predicate->first_endpoint=PARSER_REGISTRY<TV,PREDICATE<TV>>::Parse(node["first"],simulation);
    predicate->second_endpoint=PARSER_REGISTRY<TV,PREDICATE<TV>>::Parse(node["second"],simulation);
    return predicate;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template<class TV> TV RIGID_STRUCTURE_POSITION_PREDICATE<TV>::
Vector(const DATA<TV>& data,const FORCE<TV>& force)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.find("RIGID_STRUCTURE_DATA")->second);
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
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
