//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class FRAME
///////////////////////////////////////////////////////////////////////
#ifndef __FRAME__
#define __FRAME__

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
public:
    enum DEFINITIONS{STATIC_SIZE=TV::SizeAtCompileTime+Quaternion<T>::Coefficients::SizeAtCompileTime};
    TV position;
    Quaternion<T> rotation;

    FRAME(){}
    ~FRAME(){}

    const Eigen::Matrix<T,STATIC_SIZE,1> Pack()
    {
        Eigen::Matrix<T,STATIC_SIZE,1> packed;
        packed.template block<TV::SizeAtCompileTime,1>(0,0)=position;
        packed.template block<Quaternion<T>::Coefficients::SizeAtCompileTime,1>(TV::SizeAtCompileTime,0)=rotation.coeffs();
        return packed;
    }

    const void Unpack(const Eigen::Matrix<T,STATIC_SIZE,1>& packed)
    {
        std::cout<<packed<<std::endl;
        position=packed.template block<TV::SizeAtCompileTime,1>(0,0);
        rotation=packed.template block<Quaternion<T>::Coefficients::SizeAtCompileTime,1>(TV::SizeAtCompileTime,0);
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(position);
    }

};
}
#endif
