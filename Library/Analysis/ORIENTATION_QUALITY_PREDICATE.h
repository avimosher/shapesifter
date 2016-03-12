#ifndef __ORIENTATION_QUALITY_PREDICATE__
#define __ORIENTATION_QUALITY_PREDICATE__

#include <Analysis/PREDICATE.h>
namespace Mechanics{
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
