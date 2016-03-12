#ifndef __LESS_THAN_PREDICATE__
#define __LESS_THAN_PREDICATE__

#include <Analysis/PREDICATE.h>
namespace Mechanics{
template<class TV>
class LESS_THAN_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    std::shared_ptr<PREDICATE<TV>> a;
    std::shared_ptr<PREDICATE<TV>> b;

    LESS_THAN_PREDICATE() {subtype=PREDICATE<TV>::SCALAR;}

    T Scalar(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("less_than")
};
}
#endif
