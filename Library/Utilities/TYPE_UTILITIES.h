//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
#ifndef __TYPE_UTILITIES__
#define __TYPE_UTILITIES__

#include <cereal/types/polymorphic.hpp>
#include <Eigen/Dense>
using namespace Eigen;

#define GENERIC_SCALAR_TYPE_DEFINITION(TYPE,T) \
    template class TYPE<Matrix<T,1,1>>;        \
    template class TYPE<Matrix<T,1,2>>;        \
    template class TYPE<Matrix<T,1,3>>;

#define GENERIC_TYPE_DEFINITION(TYPE)          \
    GENERIC_SCALAR_TYPE_DEFINITION(TYPE,float) \
    GENERIC_SCALAR_TYPE_DEFINITION(TYPE,double)

#define GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,T) \
    CEREAL_REGISTER_TYPE(TYPE<Matrix<T,1,1> >);    \
    CEREAL_REGISTER_TYPE(TYPE<Matrix<T,1,2> >);    \
    CEREAL_REGISTER_TYPE(TYPE<Matrix<T,1,3> >);

#define GENERIC_CEREAL_REGISTRATION(TYPE) \
    GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,float) \
    GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,double)



#endif
