#ifndef __DRIVER__
#define __DRIVER__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV>
class DRIVER
{
    typedef typename TV::Scalar T;
public:
    std::shared_ptr<SIMULATION<TV>> simulation;

    DRIVER(std::shared_ptr<SIMULATION<TV>> simulation)
        :simulation(simulation)
    {}
    ~DRIVER(){}

    void Initialize();
    void Advance_One_Time_Step(const T target_time,bool &done);
    void Advance_To_Target_Time(const T target_time);
    void Execute_Main_Program();
};
}
#endif
