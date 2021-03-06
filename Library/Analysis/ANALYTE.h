#ifndef __ANALYTE__
#define __ANALYTE__

#include <Evolution/EVOLUTION_STEP.h>
namespace Mechanics{
    template<class TV> class AGGREGATOR;
    template<class TV> class PREDICATE;

template<class TV>
class ANALYTE:public EVOLUTION_STEP<TV>
{
    typedef typename TV::Scalar T;
public:
    std::string name;
    std::shared_ptr<PREDICATE<TV>> condition;
    std::shared_ptr<PREDICATE<TV>> predicate;
    std::shared_ptr<AGGREGATOR<TV>> aggregator;

    ANALYTE(){}

    void Step(SIMULATION<TV>& simulation,const T dt,const T time);
    void Finalize(Json::Value& root);
    DEFINE_TYPE_NAME("ANALYTE")
};
}
#endif
