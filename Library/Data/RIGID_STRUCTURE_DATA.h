//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE_DATA
///////////////////////////////////////////////////////////////////////
#ifndef __RIGID_STRUCTURE_DATA__
#define __RIGID_STRUCTURE_DATA__

#include <Data/DATA_TYPE.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class RIGID_STRUCTURE;

template<class TV>
class RIGID_STRUCTURE_DATA:public DATA_TYPE<TV>,public std::vector<std::unique_ptr<RIGID_STRUCTURE<TV>>>
{
    typedef typename TV::Scalar T;

public:
    RIGID_STRUCTURE_DATA();
    ~RIGID_STRUCTURE_DATA();

    DEFINE_TYPE_NAME("RIGID_STRUCTURE_DATA")
};
}
#endif
