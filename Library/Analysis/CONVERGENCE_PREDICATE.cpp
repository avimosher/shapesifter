#include <Analysis/CONVERGENCE_PREDICATE.h>
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Evolution/EVOLUTION.h>
#include <Evolution/EVOLUTION_STEP.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar CONVERGENCE_PREDICATE<TV>::
Scalar(const SIMULATION<TV>& simulation)
{
    auto solver=simulation.evolution.Find(solver_name);
    return solver->Success();
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(CONVERGENCE_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(CONVERGENCE_PREDICATE,PREDICATE)
{
    auto predicate=std::make_shared<CONVERGENCE_PREDICATE<TV>>();
    Parse_String(node["solver"],predicate->solver_name);
    return predicate;
}

