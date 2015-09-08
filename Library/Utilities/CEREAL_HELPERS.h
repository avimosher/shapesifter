#ifndef __CEREAL_HELPERS__
#define __CEREAL_HELPERS__

#include <Utilities/EIGEN_HELPERS.h>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#undef CEREAL_BIND_TO_ARCHIVES
#define CEREAL_BIND_TO_ARCHIVES(...)                                     \
    namespace cereal {                                                   \
    namespace detail {                                                   \
    template<>                                                           \
    struct init_binding<__VA_ARGS__> {                                   \
        static bind_to_archives<__VA_ARGS__> const & b;                  \
    };                                                                   \
    bind_to_archives<__VA_ARGS__> const & init_binding<__VA_ARGS__>::b = \
        ::cereal::detail::StaticObject<                                  \
            bind_to_archives<__VA_ARGS__>                                \
        >::getInstance().bind();                                         \
    }} // end namespaces

#undef CEREAL_REGISTER_TYPE
#define CEREAL_REGISTER_TYPE(...)                                 \
  namespace cereal {                                              \
  namespace detail {                                              \
  template <>                                                     \
  struct binding_name<__VA_ARGS__>                                \
  {                                                               \
    STATIC_CONSTEXPR char const * name() { return #__VA_ARGS__; } \
  };                                                              \
  } } /* end namespaces */                                        \
  CEREAL_BIND_TO_ARCHIVES(__VA_ARGS__)


namespace cereal
{
template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
typename std::enable_if<traits::is_output_serializable<BinaryData<_Scalar>, BinaryOutputArchive>::value, void>::type
save(BinaryOutputArchive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> const & m)
{
    int32_t rows = m.rows();
    int32_t cols = m.cols();
    ar(rows);
    ar(cols);
    ar(binary_data(m.data(), rows * cols * sizeof(_Scalar)));
}

template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
typename std::enable_if<traits::is_input_serializable<BinaryData<_Scalar>, BinaryInputArchive>::value, void>::type
load(BinaryInputArchive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & m)
{
    int32_t rows;
    int32_t cols;
    ar(rows);
    ar(cols);
    m.resize(rows, cols);
    ar(binary_data(m.data(), static_cast<std::size_t>(rows * cols * sizeof(_Scalar))));
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
void
save(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> const & m)
{
    int32_t rows = m.rows();
    int32_t cols = m.cols();
    ar(CEREAL_NVP(rows));
    ar(CEREAL_NVP(cols));
    std::vector<_Scalar> v;v.assign(m.data(),m.data()+rows*cols);
    ar(make_nvp("data",v));
    
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
void
load(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & m)
{
    int32_t rows;
    int32_t cols;
    ar(rows);
    ar(cols);
    m.resize(rows, cols);
    std::vector<_Scalar> v;
    ar(make_nvp("data",v));
    std::copy(v.begin(),v.end(),m.data());
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

// define cereal serializers for std::pair
/*template<class Archive,class T1,class T2>
void Serialize(Archive& archive,std::pair<T1,T2>& p)
{
    archive(p.first,p.second);
    }*/

}
#endif

