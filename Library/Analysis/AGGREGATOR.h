#ifndef __AGGREGATOR__
#define __AGGREGATOR__

#include <Utilities/TYPE_UTILITIES.h>
#include <memory>

namespace Mechanics{
template<class TV> class PREDICATE;
template<class TV> class SIMULATION;

template<class TV>
class AGGREGATOR
{
    typedef typename TV::Scalar T;
public:
    int subtype;

    AGGREGATOR(){}
    virtual ~AGGREGATOR(){}

    virtual void Aggregate_Subtype(const T aggregate){}
    virtual void Aggregate_Subtype(const TV& aggregate){}
    void Aggregate(std::shared_ptr<PREDICATE<TV>> predicate,const SIMULATION<TV>& simulation);
    virtual void Print_Report(std::ostream& out){}
};

template<class TV>
class AVERAGE_AGGREGATOR:public AGGREGATOR<TV>
{
    typedef typename TV::Scalar T;
    int t_count;
    T t_total;
    int tv_count;
    TV tv_total;
public:
    AVERAGE_AGGREGATOR();
    ~AVERAGE_AGGREGATOR(){};

    virtual void Aggregate_Subtype(const T aggregate)
    {t_count++;t_total+=aggregate;}

    virtual void Aggregate_Subtype(const TV& aggregate)
    {tv_count++;tv_total+=aggregate;}

    virtual void Print_Report(std::ostream& out);
    DEFINE_TYPE_NAME("average")
};

template<class TV>
class SUM_AGGREGATOR:public AGGREGATOR<TV>
{
    typedef typename TV::Scalar T;
    int t_count;
    T t_total;
    int tv_count;
    TV tv_total;
public:
    using AGGREGATOR<TV>::subtype;
    SUM_AGGREGATOR();

    virtual void Aggregate_Subtype(const T aggregate)
    {t_count++;t_total+=aggregate;}

    virtual void Aggregate_Subtype(const TV& aggregate)
    {tv_count++;tv_total+=aggregate;}

    virtual void Print_Report(std::ostream& out);
    DEFINE_TYPE_NAME("sum")
};

template<class TV>
class HISTOGRAM_AGGREGATOR:public AGGREGATOR<TV>
{
    typedef typename TV::Scalar T;
public:
    int number_of_bins;
    T minimum_value;
    T maximum_value;
    T bin_width;
    std::vector<int> bins;
    HISTOGRAM_AGGREGATOR();

    virtual void Aggregate_Subtype(const T aggregate)
    {int bin_index=clamp((int)((aggregate-minimum_value)/bin_width),0,number_of_bins-1);
    bins[bin_index]++;}

    virtual void Aggregate_Subtype(const TV& aggregate)
    {}

    virtual void Print_Report(std::ostream& out);
    DEFINE_TYPE_NAME("histogram")
};

template<class TV>
class RECORD_AGGREGATOR:public AGGREGATOR<TV>
{
    typedef typename TV::Scalar T;
public:
    std::vector<T> scalar_record;
    std::vector<TV> vector_record;
    RECORD_AGGREGATOR();

    virtual void Aggregate_Subtype(const T aggregate)
    {scalar_record.push_back(aggregate);}

    virtual void Aggregate_Subtype(const TV& aggregate)
    {vector_record.push_back(aggregate);}

    virtual void Print_Report(std::ostream& out);
    DEFINE_TYPE_NAME("record")
};

}
#endif
