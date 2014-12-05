//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE
///////////////////////////////////////////////////////////////////////
#ifndef __RIGID_STRUCTURE__
#define __RIGID_STRUCTURE__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{

template<class TV>
class RIGID_STRUCTURE
{
    typedef typename TV::Scalar T;

public:
#if 0
    FRAME<TV> frame;
    TWIST<TV> twist;
    MOMENT<TV> moi;
#endif

    RIGID_STRUCTURE();
    ~RIGID_STRUCTURE();
};
}
#endif
