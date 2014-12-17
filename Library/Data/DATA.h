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
#include <Eigen/Geometry>
#include <osg/Group>
#include <unordered_map>

namespace Mechanics{
class QUALITY;

template<class TV>
class DATA : public std::unordered_map<std::string,std::shared_ptr<DATA_TYPE<TV>>>
{
    typedef typename TV::Scalar T;

public:
    DATA();
    ~DATA();

    void Variables(Matrix<T,Dynamic,1>& variables);
    void Step(QUALITY& step_quality,Matrix<T,Dynamic,1> solve_result);
    T Print_All();
    void Viewer(osg::Group*& root);
};
}
#endif
