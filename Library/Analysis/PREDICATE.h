#ifndef __PREDICATE__
#define __PREDICATE__

#include <Data/RIGID_STRUCTURE.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/Dense>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV>
class PREDICATE
{
    typedef typename TV::Scalar T;
public:
    enum SUBTYPE{SCALAR,VECTOR};
    SUBTYPE subtype;

    virtual ~PREDICATE(){}
    SUBTYPE Get_Subtype()
        {return subtype;}

    virtual T Scalar(const SIMULATION<TV>& simulation){return 0;}
    virtual TV Vector(const SIMULATION<TV>& simulation){return TV();}
};
}
#endif
