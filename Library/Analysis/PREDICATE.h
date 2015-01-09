#ifndef __PREDICATE__
#define __PREDICATE__

#include <Eigen/Dense>

namespace Mechanics{
    template<class TV> class DATA;
    template<class TV> class FORCE;
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

    virtual T Scalar(const DATA<TV>& data,const FORCE<TV>& force){}
    virtual TV Vector(const DATA<TV>& data,const FORCE<TV>& force){}
};

}

#endif
