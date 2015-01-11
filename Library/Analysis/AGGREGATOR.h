#ifndef __AGGREGATOR__
#define __AGGREGATOR__

#include <Utilities/TYPE_UTILITIES.h>
#include <memory>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE;
template<class TV> class PREDICATE;

template<class TV>
class AGGREGATOR
{
    typedef typename TV::Scalar T;
public:
    AGGREGATOR(){}
    virtual ~AGGREGATOR(){}

    virtual void Aggregate_Subtype(const T aggregate){}
    virtual void Aggregate_Subtype(const TV& aggregate){}
    void Aggregate(std::shared_ptr<PREDICATE<TV>> predicate,const DATA<TV>& data,const FORCE<TV>& force);
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

}
#endif
