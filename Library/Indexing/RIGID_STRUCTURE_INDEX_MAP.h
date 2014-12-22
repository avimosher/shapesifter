#ifndef __RIGID_STRUCTURE_INDEX_MAP__
#define __RIGID_STRUCTURE_INDEX_MAP__

#include <Data/FRAME.h>

namespace Mechanics{
template<class TV>
class RIGID_STRUCTURE_INDEX_MAP
{
    typedef typename TV::Scalar T;
    enum DEFINITIONS{STATIC_SIZE=FRAME<TV>::STATIC_SIZE};
public:
    T temperature;

    RIGID_STRUCTURE_INDEX_MAP();
    ~RIGID_STRUCTURE_INDEX_MAP();


};
}
#endif
