#ifndef __PREDICATE__
#define __PREDICATE__

#include <Utilities/TYPE_UTILITIES.h>
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
    
template<class TV>
class DISTANCE_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    std::shared_ptr<PREDICATE<TV>> first_endpoint;
    std::shared_ptr<PREDICATE<TV>> second_endpoint;

    DISTANCE_PREDICATE() {subtype=PREDICATE<TV>::SCALAR;}
    ~DISTANCE_PREDICATE() {}

    T Scalar(const DATA<TV>& data,const FORCE<TV>& force);
    DEFINE_TYPE_NAME("distance")
};

template<class TV>
class RIGID_STRUCTURE_POSITION_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    std::string name;
    TV offset;

    RIGID_STRUCTURE_POSITION_PREDICATE() {subtype=PREDICATE<TV>::VECTOR;}
    ~RIGID_STRUCTURE_POSITION_PREDICATE() {}

    TV Vector(const DATA<TV>& data,const FORCE<TV>& force);
    DEFINE_TYPE_NAME("rigid_structure_position")
};

}
#endif
