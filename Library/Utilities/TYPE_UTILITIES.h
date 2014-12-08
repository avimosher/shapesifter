//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
#ifndef __TYPE_UTILITIES__
#define __TYPE_UTILITIES__

#include <cereal/types/polymorphic.hpp>
#include <Eigen/Dense>
using namespace Eigen;

#define GENERIC_SCALAR_TYPE_DEFINITION(TYPE,T) \
    template class TYPE<Matrix<T,1,1>>; \
    template class TYPE<Matrix<T,1,2>>; \
    template class TYPE<Matrix<T,1,3>>;

#define GENERIC_TYPE_DEFINITION(TYPE) \
    GENERIC_SCALAR_TYPE_DEFINITION(TYPE,float) \
    GENERIC_SCALAR_TYPE_DEFINITION(TYPE,double)

#define LOCAL_CEREAL_BIND_TO_ARCHIVES(...)                           \
    namespace cereal {                                       \
    namespace detail {                                       \
    template<>                                               \
    struct init_binding<__VA_ARGS__> {                                 \
        static bind_to_archives<__VA_ARGS__> const & b;                \
    };                                                       \
    bind_to_archives<__VA_ARGS__> const & init_binding<__VA_ARGS__>::b =         \
        ::cereal::detail::StaticObject<                      \
            bind_to_archives<__VA_ARGS__>                              \
        >::getInstance().bind();                             \
    }} // end namespaces


#define LOCAL_CEREAL_REGISTER_TYPE(...)                         \
  namespace cereal {                                    \
  namespace detail {                                    \
  template <>                                           \
  struct binding_name<__VA_ARGS__>                                \
  {                                                     \
    STATIC_CONSTEXPR char const * name() { return #__VA_ARGS__; } \
  };                                                    \
  } } /* end namespaces */                              \
  LOCAL_CEREAL_BIND_TO_ARCHIVES(__VA_ARGS__)

#define GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,T) \
    LOCAL_CEREAL_REGISTER_TYPE(TYPE<Matrix<T,1,1> >);        \
    LOCAL_CEREAL_REGISTER_TYPE(TYPE<Matrix<T,1,2> >);        \
    LOCAL_CEREAL_REGISTER_TYPE(TYPE<Matrix<T,1,3> >);

#define GENERIC_CEREAL_REGISTRATION(TYPE) \
    GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,float) \
    GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,double)



#endif
