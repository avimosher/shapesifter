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

    const Eigen::Matrix<T,STATIC_SIZE,1> Pack()
    {
        Eigen::Matrix<T,STATIC_SIZE,1> packed;
        packed.template block<TV::SizeAtCompileTime,1>(0,0)=linear;
        packed.template block<T_SPIN::SizeAtCompileTime,1>(TV::SizeAtCompileTime,0)=angular;
        return packed;
    }

    const void Unpack(const Eigen::Matrix<T,STATIC_SIZE,1>& packed)
    {
        linear=packed.template block<TV::SizeAtCompileTime,1>(0,0);
        angular=packed.template block<T_SPIN::SizeAtCompileTime,1>(TV::SizeAtCompileTime,0);
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(linear,angular);
    }
};
}
#endif
