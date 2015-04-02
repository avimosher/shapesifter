#include <Analysis/AGGREGATOR.h>
#include <Analysis/ANALYTE.h>
#include <Analysis/PREDICATE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/LOG.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void ANALYTE<TV>::
Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    if(condition && !condition->Scalar(simulation)){return;}
    aggregator->Aggregate(predicate,simulation);
    aggregator->Print_Report(LOG::cout);
    LOG::cout<<std::endl;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(ANALYTE)
DEFINE_AND_REGISTER_PARSER(ANALYTE,void)
{
    auto analyte=std::make_shared<ANALYTE<TV>>();
    analyte->aggregator=PARSER_REGISTRY<TV,AGGREGATOR<TV>>::Parse(node["aggregator"],simulation);
    if(node.isMember("condition")){analyte->condition=PARSER_REGISTRY<TV,PREDICATE<TV>>::Parse(node["condition"],simulation);}
    analyte->predicate=PARSER_REGISTRY<TV,PREDICATE<TV>>::Parse(node["predicate"],simulation);
    simulation.evolution.push_back(analyte);
    return 0;
}
