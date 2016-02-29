#include <Utilities/RANDOM.h>
#include <Utilities/TYPE_UTILITIES.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class T> std::random_device& RANDOM<T>::
Random_Device() const
{
    static std::random_device device;
    return device;
}
///////////////////////////////////////////////////////////////////////
template<class T> std::mt19937& RANDOM<T>::
Generator() const
{
    static std::mt19937 generator(Random_Device()());
    return generator;
}
///////////////////////////////////////////////////////////////////////
template<class T> std::normal_distribution<>& RANDOM<T>::
Normal_Distribution()
{
    static std::normal_distribution<> normal_distribution(0,1);
    return normal_distribution;
}
///////////////////////////////////////////////////////////////////////
template<class T> std::uniform_real_distribution<T>& RANDOM<T>::
Uniform_Distribution()
{
    static std::uniform_real_distribution<T> uniform_distribution(0,1);
    return uniform_distribution;
}
///////////////////////////////////////////////////////////////////////
template<class T> T RANDOM<T>::
Uniform(const T& a,const T& b)
{
    return Uniform_Distribution()(Generator())*(b-a)+a;
}
///////////////////////////////////////////////////////////////////////
template<class T> T RANDOM<T>::
Gaussian()
{
    return Normal_Distribution()(Generator());
}
///////////////////////////////////////////////////////////////////////
template<class T> T RANDOM<T>::
Gaussian(T mean,T stddev)
{
    return mean+stddev*Gaussian();
}
///////////////////////////////////////////////////////////////////////
namespace Mechanics{
  template class RANDOM<double>;
}
