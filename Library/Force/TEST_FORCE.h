#ifndef __TEST_FORCE__
#define __TEST_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class TEST_FORCE : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
public:
    TEST_FORCE(){}
    ~TEST_FORCE(){}

    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    DEFINE_TYPE_NAME("TEST_FORCE")
};
}
#endif
