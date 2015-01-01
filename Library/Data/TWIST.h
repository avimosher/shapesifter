#ifndef __TWIST__
#define __TWIST__

#include <Data/ROTATION.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <cereal/archives/binary.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace Mechanics{

template<class TV>
class TWIST
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{STATIC_SIZE=TV::SizeAtCompileTime+T_SPIN::SizeAtCompileTime};
    TV linear;
    T_SPIN angular;

    TWIST(){}
    ~TWIST(){}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(linear,angular);
    }

};
}
#endif
