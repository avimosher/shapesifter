#ifndef __RIGID_STRUCTURE_POSITION_PREDICATE__
#define __RIGID_STRUCTURE_POSITION_PREDICATE__

#include <Analysis/PREDICATE.h>
namespace Mechanics{
template<class TV>
class RIGID_STRUCTURE_POSITION_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    std::string name;
    TV offset;

    RIGID_STRUCTURE_POSITION_PREDICATE() {subtype=PREDICATE<TV>::VECTOR;}

    TV Vector(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("rigid_structure_position")
};
}
#endif
