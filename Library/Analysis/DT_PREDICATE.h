#ifndef __DT_PREDICATE__
#define __DT_PREDICATE__

#include <Analysis/PREDICATE.h>
namespace Mechanics{
template<class TV>
class DT_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;

    DT_PREDICATE() {subtype=PREDICATE<TV>::SCALAR;}

    T Scalar(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("dt")
};
}
#endif
