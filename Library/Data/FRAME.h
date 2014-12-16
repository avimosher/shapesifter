//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class FRAME
///////////////////////////////////////////////////////////////////////
#ifndef __FRAME__
#define __FRAME__

#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/Geometry>

namespace Mechanics{

template<class TV>
class FRAME
{
    typedef typename TV::Scalar T;
public:
    TV position;
    Quaternion<T> orientation;

    FRAME();
    ~FRAME();
};
}
#endif
