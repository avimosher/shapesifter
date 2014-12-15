///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class DATA
///////////////////////////////////////////////////////////////////////
#ifndef __DATA__
#define __DATA__

#include <Data/DATA_TYPE.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <memory>
#include <vector>
#include <osg/Group>

namespace Mechanics{
class QUALITY;

template<class TV>
class DATA : public std::vector<std::unique_ptr<DATA_TYPE<TV>>>
{
    typedef typename TV::Scalar T;

public:
    bool restart;
    std::string output_directory;

    DATA();
    ~DATA();

    void Variables(Matrix<T,Dynamic,1>& variables);
    void Step(QUALITY& step_quality,Matrix<T,Dynamic,1> solve_result);
    void Write(const int frame);
    void Read(const int frame);
    T Print_All();
    void Viewer(osg::Group*& root);
};
}
#endif
