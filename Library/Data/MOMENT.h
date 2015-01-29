#ifndef __MOMENT__
#define __MOMENT__

#include <Data/ROTATION.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <cereal/archives/binary.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace Mechanics{

template<class TV>
class MOMENT
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    TV translation;
    T_SPIN rotation;

    MOMENT(){}
    ~MOMENT(){}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(translation,rotation);
    }
};
}
#endif
