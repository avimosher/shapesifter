#ifndef __QUALITY__
#define __QUALITY__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{

template<class T>
class QUALITY
{
public:
    T step_scaling;
    T last_norm;

    QUALITY()
        :step_scaling((T).25),last_norm(FLT_MAX)
    {};
    ~QUALITY(){};

    void Update(T norm){
        T quality_ratio=(last_norm-norm)/last_norm;
        if(quality_ratio>.75*step_scaling){step_scaling=std::min(.75,1.2*step_scaling);}
        else if(quality_ratio<.3*step_scaling){step_scaling=std::max(.1,.8*step_scaling);}
        last_norm=norm;
    }

    T Scale_Result(){return step_scaling;}
};
}
#endif
