#include <Analysis/DISTANCE_PREDICATE.h>
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar DISTANCE_PREDICATE<TV>::
Scalar(const SIMULATION<TV>& simulation)
{
    return simulation.data.Minimum_Offset(first_endpoint->Vector(simulation),second_endpoint->Vector(simulation)).norm();
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

