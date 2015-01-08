//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE
///////////////////////////////////////////////////////////////////////
#ifndef __RIGID_STRUCTURE__
#define __RIGID_STRUCTURE__

#include <Data/FRAME.h>
#include <Data/TWIST.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <cereal/archives/binary.hpp>

namespace Mechanics{

template<class TV>
class RIGID_STRUCTURE
{
    typedef typename TV::Scalar T;

public:
    enum DEFINITIONS{STATIC_SIZE=FRAME<TV>::STATIC_SIZE};
    std::string name;
    FRAME<TV> frame;
    TWIST<TV> twist;
#if 0
    MOMENT<TV> moi;
#endif

    RIGID_STRUCTURE();
    ~RIGID_STRUCTURE(){};

    template<class Archive>
    void serialize(Archive& archive) {
        archive(frame,twist);
    }

    DEFINE_TYPE_NAME("RIGID_STRUCTURE")
};
}
#endif
