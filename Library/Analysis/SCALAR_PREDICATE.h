#ifndef __SCALAR_PREDICATE__
#define __SCALAR_PREDICATE__

#include <Analysis/PREDICATE.h>
namespace Mechanics{
template<class TV>
class SCALAR_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    T scalar;

    SCALAR_PREDICATE() {subtype=PREDICATE<TV>::SCALAR;}

    T Scalar(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("scalar")
};
}
#endif
