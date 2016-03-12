#include <Analysis/DT_PREDICATE.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar DT_PREDICATE<TV>::
Scalar(const SIMULATION<TV>& simulation)
{
    return simulation.dt;
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(DT_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(DT_PREDICATE,PREDICATE)
{
    return std::make_shared<DT_PREDICATE<TV>>();
}
