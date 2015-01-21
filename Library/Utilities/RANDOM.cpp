#include <Utilities/RANDOM.h>
#include <Utilities/TYPE_UTILITIES.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> std::random_device& RANDOM<TV>::
Random_Device()
{
    static std::random_device device;
    return device;
}
///////////////////////////////////////////////////////////////////////
template<class TV> std::mt19937& RANDOM<TV>::
Generator()
{
    static std::mt19937 generator(Random_Device()());
    return generator;
}
///////////////////////////////////////////////////////////////////////
template<class TV> std::normal_distribution<>& RANDOM<TV>::
Normal_Distribution()
{
    static std::normal_distribution<> normal_distribution(0,1);
    return normal_distribution;
}
///////////////////////////////////////////////////////////////////////
template<class TV> std::uniform_real_distribution<typename TV::Scalar>& RANDOM<TV>::
Uniform_Distribution()
{
    static std::uniform_real_distribution<T> uniform_distribution(0,1);
    return uniform_distribution;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar RANDOM<TV>::
Uniform(const T& a,const T& b)
{
    return Uniform_Distribution()(Generator())*(b-a)+a;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar RANDOM<TV>::
Gaussian()
{
    return Normal_Distribution()(Generator());
}
///////////////////////////////////////////////////////////////////////
template<class TV> TV RANDOM<TV>::
Direction()
{
    TV direction;
    for(int i=0;i<TV::RowsAtCompileTime;i++){direction(i)=Gaussian();}
    return direction.normalized();
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RANDOM)
