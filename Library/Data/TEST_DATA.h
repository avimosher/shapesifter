//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class TEST_DATA
///////////////////////////////////////////////////////////////////////
#ifndef __TEST_DATA__
#define __TEST_DATA__

#include <Data/DATA_TYPE.h>
#include <iostream>
#include <cereal/archives/json.hpp>

namespace Mechanics{

template<class TV>
class TEST_DATA:public DATA_TYPE<TV>
{
    typedef typename TV::Scalar T;

    T internal_data; // v. simple
public:
    TEST_DATA(){internal_data=5;}
    ~TEST_DATA(){}

    virtual int Size(){
        return 1;
    }

    virtual Matrix<T,Dynamic,1> Variables() {
        Matrix<T,Dynamic,1> variables(1,1);
        variables(0,0)=internal_data;
        return variables;
    }

    virtual void Step(const Matrix<T,Dynamic,1>& variables) {
        std::cout<<"Stepped"<<std::endl;
        internal_data=variables(0,0);
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(cereal::make_nvp("internal_data",internal_data));
    }
};
}
#endif
