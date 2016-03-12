#ifndef __VECTOR_PREDICATE__
#define __VECTOR_PREDICATE__

#include <Analysis/PREDICATE.h>
namespace Mechanics{
template<class TV>
class VECTOR_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    TV vector;

    VECTOR_PREDICATE() {subtype=PREDICATE<TV>::VECTOR;}

    TV Vector(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("vector")
};
}
#endif
