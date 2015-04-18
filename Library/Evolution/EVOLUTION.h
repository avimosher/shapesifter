#ifndef __EVOLUTION__
#define __EVOLUTION__

#include <Utilities/TYPE_UTILITIES.h>
#include <memory>

namespace Mechanics{
template<class TV> class EVOLUTION_STEP;
template<class TV> class SIMULATION;

template<class TV>
class EVOLUTION : public std::vector<std::shared_ptr<EVOLUTION_STEP<TV>>>
{
    typedef typename TV::Scalar T;

public:
    EVOLUTION(){}
    ~EVOLUTION(){}

    std::shared_ptr<EVOLUTION_STEP<TV>> Find(const std::string& name) const {
        Finder<std::shared_ptr<EVOLUTION_STEP<TV>>> finder={name};
        auto found=std::find_if(this->begin(),this->end(),finder);
        if(found==this->end()){std::cout<<"Could not find evolution step: "<<name<<std::endl;exit(-1);}
        return *found;
    }

    T Compute_Dt(SIMULATION<TV>& simulation,const T time,const T target_time,bool& done);
    void Advance_One_Step(SIMULATION<TV>& simulation,const T dt,const T time);
    void Finalize();
};
}
#endif
