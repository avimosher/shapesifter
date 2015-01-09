#include <Analysis/ANALYTE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void ANALYTE<TV>::
Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(ANALYTE)
DEFINE_AND_REGISTER_PARSER(ANALYTE,void)
{
    auto analyte=std::make_shared<ANALYTE<TV>>();
    //analyte->aggregator=AGGREGATOR_PARSER_REGISTRY<TV>::Parse(node["aggregator"]);
    //analyte->condition=CONDITION_PARSER_REGISTRY<TV>::Parse(node["condition"]);
    //analyte->predicate=PREDICATE_PARSER_REGISTRY<TV>::Parse(node["predicate"]);
    simulation.evolution.push_back(analyte);
    return 0;
}
