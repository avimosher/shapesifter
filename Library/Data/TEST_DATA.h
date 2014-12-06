//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class TEST_DATA
///////////////////////////////////////////////////////////////////////
#ifndef __TEST_DATA__
#define __TEST_DATA__

#include <Data/DATA_TYPE.h>

namespace Mechanics{

template<class TV>
class TEST_DATA:public DATA_TYPE<TV>
{
    typedef typename TV::Scalar T;

    T internal_data; // v. simple
public:
    TEST_DATA(){internal_data=5;}
    ~TEST_DATA(){}

    virtual Matrix<T,Dynamic,1> Variables() {
        Matrix<T,Dynamic,1> variables(1,1);
        variables(0,0)=internal_data;
        return variables;
    };
};
}
#endif
