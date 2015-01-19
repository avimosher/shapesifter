#ifndef __RANDOM__
#define __RANDOM__

#include <random>

namespace Mechanics{

template<class TV>
class RANDOM
{
    typedef typename TV::Scalar T;

public:
    std::random_device& Random_Device();
    std::mt19937& Generator();
    std::normal_distribution<>& Normal_Distribution();
    T Gaussian();
    TV Direction();
};
}
#endif
