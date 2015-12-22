#ifndef __RANDOM__
#define __RANDOM__

#include <Utilities/CEREAL_HELPERS.h>
#include <random>

namespace Mechanics{

template<class T>
class RANDOM
{
public:
    std::random_device& Random_Device() const;
    std::mt19937& Generator() const;
    std::normal_distribution<>& Normal_Distribution();
    std::uniform_real_distribution<T>& Uniform_Distribution();
    T Uniform(const T& a,const T& b);
    T Gaussian();
    T Gaussian(T mean,T stddev);

    template<class Archive>
    void save(Archive& archive) const
    {
        std::stringstream mt_stringstream;mt_stringstream<<Generator();
        archive(mt_stringstream.str());
    }

    template<class Archive>
    void load(Archive& archive)
    {
        std::string mt_string;
        archive(mt_string);
        std::stringstream mt_stringstream(mt_string);
        mt_stringstream>>Generator();
    }

    template<class TV> TV Direction(){
        TV direction;
        for(int i=0;i<TV::RowsAtCompileTime;i++){direction(i)=Gaussian();}
        return direction.normalized();
    }

    template<class TV> void Direction(TV& direction){
        for(int i=0;i<direction.size();i++){direction(i)=Gaussian();}
        direction.normalize();
    }
};
}
#endif
