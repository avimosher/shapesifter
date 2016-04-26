#ifndef __CONVERGENCE_PREDICATE__
#define __CONVERGENCE_PREDICATE__

#include <Analysis/PREDICATE.h>
namespace Mechanics{
template<class TV>
class CONVERGENCE_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    std::string solver_name;

    CONVERGENCE_PREDICATE() {subtype=PREDICATE<TV>::SCALAR;}

    T Scalar(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("convergence")
};
}
#endif
