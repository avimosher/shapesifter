#include <Analysis/AGGREGATOR.h>
#include <Analysis/ANALYTE.h>
#include <Analysis/PREDICATE.h>
#include <Driver/SIMULATION.h>
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
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ANALYTE<TV>::
Finalize()
{
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
    analyte->aggregator->subtype=analyte->predicate->subtype;
    Json::Value prerequisites=node["prerequisites"];
    if(!prerequisites.isNull()){
        for(Json::ValueIterator it=prerequisites.begin();it!=prerequisites.end();it++){
            analyte->prerequisites.push_back((*it).asString());}}
    simulation.evolution.push_back(analyte);
    return 0;
}
