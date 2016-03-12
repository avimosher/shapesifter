#include <Analysis/VECTOR_PREDICATE.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> TV VECTOR_PREDICATE<TV>::
Vector(const SIMULATION<TV>& simulation)
{
    return vector;
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(VECTOR_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(VECTOR_PREDICATE,PREDICATE)
{
    auto predicate=std::make_shared<VECTOR_PREDICATE<TV>>();
    Parse_Vector(node["vector"],predicate->vector);
    return predicate;
}
