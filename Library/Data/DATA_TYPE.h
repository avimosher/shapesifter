//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class DATA_TYPE
///////////////////////////////////////////////////////////////////////
#ifndef __DATA_TYPE__
#define __DATA_TYPE__

#include <Data/DATA_TYPE.h>
#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{

template<class TV>
class DATA_TYPE
{
    typedef typename TV::Scalar T;

public:
    DATA_TYPE();
    ~DATA_TYPE();
    virtual Matrix<T,Dynamic,1> Variables()=0;
};
}
#endif
