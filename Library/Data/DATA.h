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
    DATA(){}
    ~DATA(){}

    TV Wrap(const TV& unwrapped) const{
        return unwrapped; // TODO: wrap to domain boundaries
    }

    TV Minimum_Offset(const TV& X1,const TV& X2) const{
        return X2-X1;
    }

    int Velocity_DOF() const;
    void Pack_Velocities(Matrix<T,Dynamic,1>& velocities);
    void Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities);
    void Pack_Positions(Matrix<T,Dynamic,1>& positions);
    void Unpack_Positions(const Matrix<T,Dynamic,1>& positions);
    void Step(QUALITY& step_quality);
    void Viewer(osg::Group*& root);
};
}
#endif
