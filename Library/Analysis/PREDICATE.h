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

    virtual T Scalar(const SIMULATION<TV>& simulation){}
    virtual TV Vector(const SIMULATION<TV>& simulation){}
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

    T Scalar(const SIMULATION<TV>& simulation);
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

    TV Vector(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("rigid_structure_position")
};

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

template<class TV>
class ORIENTATION_QUALITY_PREDICATE:public PREDICATE<TV>
{
    typedef typename TV::Scalar T;
public:
    using PREDICATE<TV>::subtype;
    RIGID_STRUCTURE<TV> receptor;
    std::string binder_name;
    std::string occluder_name;
    T distance_limit;
    T around_bond_angle_limit;
    T out_of_bond_angle_limit;
    TV target_bond_orientation;
    TV binder_bond_vector;
    T target_height;
    
    ORIENTATION_QUALITY_PREDICATE() {subtype=PREDICATE<TV>::SCALAR;}

    T Sample_Point(const TV& point,const TV& offset,const DATA<TV>& data,RIGID_STRUCTURE<TV>& receptor,const RIGID_STRUCTURE<TV>& occluder,const T height,const T one_over_distance_limit_squared);
    T Scalar(const SIMULATION<TV>& simulation);
    DEFINE_TYPE_NAME("orientation_quality");
};

}
#endif
