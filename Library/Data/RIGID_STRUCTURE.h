//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE
///////////////////////////////////////////////////////////////////////
#ifndef __RIGID_STRUCTURE__
#define __RIGID_STRUCTURE__

namespace Mechanics{

template<class TV>
class RIGID_STRUCTURE
{
    typedef typename TV::Scalar T;

public:
    FRAME<TV> frame;
    TWIST<TV> twist;
    MOMENT<TV> moi;

    RIGID_STRUCTURE();
    ~RIGID_STRUCTURE();
};
}
#endif
