#ifndef __DISTANCE_PREDICATE__
#define __DISTANCE_PREDICATE__

#include <Analysis/PREDICATE.h>
namespace Mechanics{
template<class TV>
class DISTANCE_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    std::shared_ptr<PREDICATE<TV>> first_endpoint;
    std::shared_ptr<PREDICATE<TV>> second_endpoint;

    DISTANCE_PREDICATE() {subtype=PREDICATE<TV>::SCALAR;}

    T Scalar(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("distance")
};
}
#endif
