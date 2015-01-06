#ifndef __CEREAL_HELPERS__
#define __CEREAL_HELPERS__

#include <Utilities/EIGEN_HELPERS.h>
#include <cereal/archives/binary.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

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

// define cereal serializers for Rotation1D, Rotation2D and Quaternion
  template<class Archive,class T>
      void serialize(Archive& archive,Eigen::Rotation1D<T>& rotation)
  {}
  
  template<class Archive,class T>
      void serialize(Archive& archive,Eigen::Rotation2D<T>& rotation)
  {
      archive(rotation.angle());
  }

  template<class Archive,class T>
      void serialize(Archive& archive,Eigen::Quaternion<T>& rotation)
  {
      archive(rotation.coeffs());
  }

}


#endif

