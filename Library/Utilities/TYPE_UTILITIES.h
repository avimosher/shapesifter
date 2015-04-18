#ifndef __TYPE_UTILITIES__
#define __TYPE_UTILITIES__

#include <Utilities/CEREAL_HELPERS.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;

#define GENERIC_SCALAR_TYPE_DEFINITION(TYPE,T) \
    template class TYPE<Matrix<T,3,1>>;

#define GENERIC_TYPE_DEFINITION(TYPE)          \
    GENERIC_SCALAR_TYPE_DEFINITION(TYPE,float) \
    GENERIC_SCALAR_TYPE_DEFINITION(TYPE,double)

#define GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,T) \
    CEREAL_REGISTER_TYPE(TYPE<Matrix<T,3,1> >);

#define GENERIC_CEREAL_REGISTRATION(TYPE) \
    GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,float) \
    GENERIC_SCALAR_CEREAL_REGISTRATION(TYPE,double)

#define DEFINE_TYPE_NAME(name)       \
    static std::string Static_Name() \
    {return name;}                   \
                                     \
    virtual std::string Name()       \
    {return Static_Name();}

template<class T>
struct Finder
{
    const std::string name;
    
    bool operator()(const T& other) const
    {
        return other->Name()==name;
    }
};

#endif
