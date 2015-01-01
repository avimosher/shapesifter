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

    TV Wrap(const TV& unwrapped) const{
        return unwrapped; // TODO: wrap to domain boundaries
    }

    TV Minimum_Offset(const TV& X1,const TV& X2) const{
        return X2-X1;
    }

    void Variables(Matrix<T,Dynamic,1>& variables);
    void Step(QUALITY& step_quality,Matrix<T,Dynamic,1> solve_result);
    T Print_All();
    void Viewer(osg::Group*& root);
};
}
#endif
