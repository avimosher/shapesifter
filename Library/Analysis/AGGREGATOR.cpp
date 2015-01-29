#include <Analysis/AGGREGATOR.h>
#include <Analysis/PREDICATE.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template<class TV> void AGGREGATOR<TV>::
Aggregate(std::shared_ptr<PREDICATE<TV>> predicate,const SIMULATION<TV>& simulation)
{
    switch(predicate->Get_Subtype()){
        case PREDICATE<TV>::SCALAR:
            Aggregate_Subtype(predicate->Scalar(simulation));
            break;
        case PREDICATE<TV>::VECTOR:
            Aggregate_Subtype(predicate->Vector(simulation));
            break;
        default:
            // TODO: error
            break;
    }
}
GENERIC_TYPE_DEFINITION(AGGREGATOR)
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template<class TV> AVERAGE_AGGREGATOR<TV>::
AVERAGE_AGGREGATOR()
:t_count(0),t_total(0),tv_count(0),tv_total(TV())
{}
//////////////////////////////////////////////////////////////////////
template<class TV> void AVERAGE_AGGREGATOR<TV>::
Print_Report(std::ostream& out)
{
    if(t_count>0){
        out<<t_total/t_count;
    }
    if(tv_count>0){
        out<<tv_total/tv_count;
    }
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(AVERAGE_AGGREGATOR)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(AVERAGE_AGGREGATOR,AGGREGATOR)
{
    auto aggregator=std::make_shared<AVERAGE_AGGREGATOR<TV>>();
    return aggregator;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template<class TV> SUM_AGGREGATOR<TV>::
SUM_AGGREGATOR()
    :t_count(0),t_total(0),tv_count(0),tv_total(TV())
{}
//////////////////////////////////////////////////////////////////////
template<class TV> void SUM_AGGREGATOR<TV>::
Print_Report(std::ostream& out)
{
    if(t_count>0){
        out<<t_total;
    }
    if(tv_count>0){
        out<<tv_total;
    }
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SUM_AGGREGATOR)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(SUM_AGGREGATOR,AGGREGATOR)
{
    auto aggregator=std::make_shared<SUM_AGGREGATOR<TV>>();
    return aggregator;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
