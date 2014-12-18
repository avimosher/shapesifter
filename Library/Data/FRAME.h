//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class FRAME
///////////////////////////////////////////////////////////////////////
#ifndef __FRAME__
#define __FRAME__

#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <cereal/archives/binary.hpp>

namespace cereal
{
  template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
  typename std::enable_if<traits::is_output_serializable<BinaryData<_Scalar>, Archive>::value, void>::type
    save(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> const & m)
    {
      int32_t rows = m.rows();
      int32_t cols = m.cols();
      ar(rows);
      ar(cols);
      ar(binary_data(m.data(), rows * cols * sizeof(_Scalar)));
    }

  template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
  typename std::enable_if<traits::is_input_serializable<BinaryData<_Scalar>, Archive>::value, void>::type
  load(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & m)
    {
      int32_t rows;
      int32_t cols;
      ar(rows);
      ar(cols);

      m.resize(rows, cols);

      ar(binary_data(m.data(), static_cast<std::size_t>(rows * cols * sizeof(_Scalar))));
    }
}

namespace Mechanics{
template<typename T> struct incomplete;

template<class TV>
class FRAME
{
    typedef typename TV::Scalar T;
public:
    enum DEFINITIONS{STATIC_SIZE=TV::SizeAtCompileTime+Quaternion<T>::Coefficients::SizeAtCompileTime};
    TV position;
    Quaternion<T> orientation;

    FRAME(){}
    ~FRAME(){}

    const Eigen::Matrix<T,STATIC_SIZE,1> Pack()
    {
        Eigen::Matrix<T,STATIC_SIZE,1> packed;
        packed.template block<TV::SizeAtCompileTime,1>(0,0)=position;
        packed.template block<Quaternion<T>::Coefficients::SizeAtCompileTime,1>(TV::SizeAtCompileTime,0)=orientation.coeffs();
        return packed;
    }

    const void Unpack(const Eigen::Matrix<T,STATIC_SIZE,1>& packed)
    {
        std::cout<<packed<<std::endl;
        std::cout<<"Position before: "<<position<<std::endl;
        position=packed.template block<TV::SizeAtCompileTime,1>(0,0);
        orientation=packed.template block<Quaternion<T>::Coefficients::SizeAtCompileTime,1>(TV::SizeAtCompileTime,0);
        std::cout<<"Position after: "<<position<<std::endl;
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(position);
    }

};
}
#endif
