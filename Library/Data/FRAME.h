#ifndef __FRAME__
#define __FRAME__

#include <Data/ROTATION.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <cereal/archives/binary.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace Mechanics{

template<class TV>
class FRAME
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::ORIENTATION T_ORIENTATION;
public:
    enum DEFINITIONS{STATIC_SIZE=TV::SizeAtCompileTime+ROTATION<TV>::SizeAtCompileTime};
    TV position;
    T_ORIENTATION orientation;

    FRAME()
        :orientation(T_ORIENTATION::Identity())
    {}
    ~FRAME(){}

    TV operator*(const TV& vector) const
    {
        return orientation._transformVector(vector)+position;
    }

    const Eigen::Matrix<T,STATIC_SIZE,1> Pack()
    {
        Eigen::Matrix<T,STATIC_SIZE,1> packed;
        packed.template block<TV::SizeAtCompileTime,1>(0,0)=position;
        packed.template block<ROTATION<TV>::SizeAtCompileTime,1>(TV::SizeAtCompileTime,0)=Get_Rotation_Coefficients(orientation);
        return packed;
    }

    const void Unpack(const Eigen::Matrix<T,STATIC_SIZE,1>& packed)
    {
        position=packed.template block<TV::SizeAtCompileTime,1>(0,0);
        Set_Rotation_Coefficients(orientation,packed.template block<ROTATION<TV>::SizeAtCompileTime,1>(TV::SizeAtCompileTime,0));
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(position,orientation);
    }

};
}
#endif
