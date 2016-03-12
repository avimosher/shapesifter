#include <Analysis/SCALAR_PREDICATE.h>
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar SCALAR_PREDICATE<TV>::
Scalar(const SIMULATION<TV>& simulation)
{
    return scalar;
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SCALAR_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(SCALAR_PREDICATE,PREDICATE)
{
    auto predicate=std::make_shared<SCALAR_PREDICATE<TV>>();
    Parse_Scalar(node["scalar"],predicate->scalar);
    return predicate;
}

