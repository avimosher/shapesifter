#ifndef __DEFINED_FORCE__
#define __DEFINED_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class DEFINED_FORCE:public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    struct FORCE_POKE{
        int s;
        TV v;
        TV f;
    };
    std::vector<FORCE_POKE> forces;

    DEFINED_FORCE(){}
    ~DEFINED_FORCE(){}

    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    DEFINE_TYPE_NAME("DEFINED_FORCE")
};
}
#endif
