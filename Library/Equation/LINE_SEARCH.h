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

    LINE_SEARCH();
    ~LINE_SEARCH(){}

    void Step(SIMULATION<TV>& simulation,const T dt,const T time);

};
}
#endif
