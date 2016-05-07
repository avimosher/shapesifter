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
            LOG::cout<<"Attempting to aggregate invalid predicate type"<<std::endl;
            exit(EXIT_FAILURE);
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
Print_Report(Json::Value& node)
{
    if(subtype==PREDICATE<TV>::SCALAR){
        node=Json::Value(t_total/std::max(1,t_count));}
    else{Set_Vector(tv_total/std::max(1,tv_count),node);}
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
Print_Report(Json::Value& node)
{
    if(subtype==PREDICATE<TV>::SCALAR){
        node=Json::Value(t_total);}
    else{Set_Vector(tv_total,node);}
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
template<class TV> HISTOGRAM_AGGREGATOR<TV>::
HISTOGRAM_AGGREGATOR()
{}
//////////////////////////////////////////////////////////////////////
template<class TV> void HISTOGRAM_AGGREGATOR<TV>::
Print_Report(Json::Value& node)
{
    node=Json::Value(Json::arrayValue);
    //std::cout<<bins[1];
    for(int i=0;i<number_of_bins;i++){node.append(bins[i]);}
    //std::cout<<"\t"<<bins[i]<<std::endl;
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(HISTOGRAM_AGGREGATOR)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(HISTOGRAM_AGGREGATOR,AGGREGATOR)
{
    auto aggregator=std::make_shared<HISTOGRAM_AGGREGATOR<TV>>();
    Parse_Scalar(node["bins"],aggregator->number_of_bins);
    Parse_Scalar(node["min"],aggregator->minimum_value);
    Parse_Scalar(node["max"],aggregator->maximum_value);
    aggregator->bins.resize(aggregator->number_of_bins);
    aggregator->bin_width=(aggregator->maximum_value-aggregator->minimum_value)/aggregator->number_of_bins;
    return aggregator;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template<class TV> RECORD_AGGREGATOR<TV>::
RECORD_AGGREGATOR()
{}
//////////////////////////////////////////////////////////////////////
template<class TV> void RECORD_AGGREGATOR<TV>::
Print_Report(Json::Value& node)
{
    node=Json::Value(Json::arrayValue);
    if(subtype==PREDICATE<TV>::SCALAR){
        for(int i=0;i<scalar_record.size();i++){node.append(scalar_record[i]);}}
    else{
        for(int i=0;i<vector_record.size();i++){
            Json::Value element;
            Set_Vector(vector_record[i],element);
            node.append(element);}}
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RECORD_AGGREGATOR)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(RECORD_AGGREGATOR,AGGREGATOR)
{
    auto aggregator=std::make_shared<RECORD_AGGREGATOR<TV>>();
    return aggregator;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

