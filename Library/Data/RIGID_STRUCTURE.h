//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE
///////////////////////////////////////////////////////////////////////
#ifndef __RIGID_STRUCTURE__
#define __RIGID_STRUCTURE__

#include <Data/FRAME.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <cereal/archives/binary.hpp>

namespace Mechanics{

template<class TV>
class RIGID_STRUCTURE
{
    typedef typename TV::Scalar T;

public:
    enum DEFINITIONS{STATIC_SIZE=FRAME<TV>::STATIC_SIZE};
    FRAME<TV> frame;
#if 0
    TWIST<TV> twist;
    MOMENT<TV> moi;
#endif

    RIGID_STRUCTURE();
    ~RIGID_STRUCTURE(){};

    template<class Archive>
    void serialize(Archive& archive) {
        std::cout<<"Rigid serialize"<<std::endl;
        archive(frame);
        std::cout<<frame.position<<std::endl;
    }

    DEFINE_TYPE_NAME("RIGID_STRUCTURE")
};
}
#endif
