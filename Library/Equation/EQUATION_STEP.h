#ifndef __EQUATION_STEP__
#define __EQUATION_STEP__

#include <Evolution/EVOLUTION_STEP.h>

namespace Mechanics{
    template<class TV> class EQUATION;
    
template<class TV>
class EQUATION_STEP:public EVOLUTION_STEP<TV>
{
    typedef typename TV::Scalar T;

public:
    EQUATION<TV>* equation;

    EQUATION_STEP(){}
    ~EQUATION_STEP(){}

    void Step(SIMULATION<TV>& simulation,const T dt,const T time);
};
}
#endif
