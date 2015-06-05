#ifndef __LINE_SEARCH__
#define __LINE_SEARCH__

#include <Evolution/EVOLUTION_STEP.h>

namespace Mechanics{
template<class TV> class EQUATION;

template<class TV>
class LINE_SEARCH:public EVOLUTION_STEP<TV>
{
    typedef typename TV::Scalar T;
public:
    std::shared_ptr<EQUATION<TV>> equation;
    LINE_SEARCH();
    ~LINE_SEARCH(){}

    void Step(SIMULATION<TV>& simulation,const T dt,const T time);
    DEFINE_TYPE_NAME("LINE_SEARCH")
};
}
#endif
