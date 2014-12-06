///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class DATA
///////////////////////////////////////////////////////////////////////
#ifndef __DATA__
#define __DATA__

#include <Data/DATA_TYPE.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <vector>

namespace Mechanics{
class QUALITY;

template<class TV>
class DATA : public std::vector<DATA_TYPE<TV>*>
{
    typedef typename TV::Scalar T;

public:
    bool restart;

    DATA();
    ~DATA();

    void Variables(Matrix<T,Dynamic,1>& variables);
    void Step(QUALITY& step_quality,Matrix<T,Dynamic,1> solve_result);
};
}
#endif
