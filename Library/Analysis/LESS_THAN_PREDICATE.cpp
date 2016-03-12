#include <Analysis/LESS_THAN_PREDICATE.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar LESS_THAN_PREDICATE<TV>::
Scalar(const SIMULATION<TV>& simulation)
{
    return a->Scalar(simulation)<b->Scalar(simulation);
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(LESS_THAN_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(LESS_THAN_PREDICATE,PREDICATE)
{
    auto predicate=std::make_shared<LESS_THAN_PREDICATE<TV>>();
    predicate->a=PARSER_REGISTRY<TV,PREDICATE<TV>>::Parse(node["a"],simulation);
    predicate->b=PARSER_REGISTRY<TV,PREDICATE<TV>>::Parse(node["b"],simulation);
    return predicate;
}
