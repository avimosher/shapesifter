#include <Analysis/AGGREGATOR.h>
#include <Analysis/ANALYTE.h>
#include <Analysis/PREDICATE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void ANALYTE<TV>::
Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    if(condition && !condition->Scalar(data,force)){return;}
    aggregator->Aggregate(predicate,data,force);
    aggregator->Print_Report(std::cout);
    std::cout<<std::endl;
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
