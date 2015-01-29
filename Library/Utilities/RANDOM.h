#ifndef __RANDOM__
#define __RANDOM__

#include <random>

namespace Mechanics{

template<class T>
class RANDOM
{
public:
    std::random_device& Random_Device();
    std::mt19937& Generator();
    std::normal_distribution<>& Normal_Distribution();
    std::uniform_real_distribution<T>& Uniform_Distribution();
    T Uniform(const T& a,const T& b);
    T Gaussian();
    T Gaussian(T mean,T stddev);

    template<class TV> TV Direction(){
        TV direction;
        for(int i=0;i<TV::RowsAtCompileTime;i++){direction(i)=Gaussian();}
        return direction.normalized();
    }
};
}
#endif
