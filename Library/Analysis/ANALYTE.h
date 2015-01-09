#ifndef __ANALYTE__
#define __ANALYTE__

#include <Evolution/EVOLUTION_STEP.h>
namespace Mechanics{
    template<class TV> class AGGREGATOR;
    template<class TV> class CONDITION;
    template<class TV> class PREDICATE;

template<class TV>
class ANALYTE:public EVOLUTION_STEP<TV>
{
    typedef typename TV::Scalar T;
public:
    CONDITION<TV>* condition;
    PREDICATE<TV>* predicate;
    AGGREGATOR<TV>* aggregator;

    ANALYTE(){}

    void Step(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time);
    DEFINE_TYPE_NAME("ANALYTE")
};
}
#endif
