#ifndef __DATA__
#define __DATA__

#include <Data/DATA_TYPE.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <memory>
#include <Eigen/Geometry>
#include <osg/Group>
#include <unordered_map>

namespace Mechanics{

template<class TV>
class DATA:public std::vector<std::shared_ptr<DATA_TYPE<TV>>>
{
    typedef typename TV::Scalar T;

public:
    std::unordered_map<std::string,T> globals;

    DATA(){}
    ~DATA(){}

    TV Wrap(const TV& unwrapped) const{
        return unwrapped; // TODO: wrap to domain boundaries
    }

    TV Minimum_Offset(const TV& X1,const TV& X2) const{
        return X2-X1;
    }

    std::shared_ptr<DATA_TYPE<TV>> Find(const std::string& name) const{
        Finder<std::shared_ptr<DATA_TYPE<TV>>> finder={name};
        return *(std::find_if(this->begin(),this->end(),finder));
    }

    int Velocity_DOF() const;
    void Pack_Velocities(Matrix<T,Dynamic,1>& velocities);
    void Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities);
    void Pack_Positions(Matrix<T,Dynamic,1>& positions);
    void Unpack_Positions(const Matrix<T,Dynamic,1>& positions);
    void Step();
    void Viewer(osg::Group*& root);

    DEFINE_TYPE_NAME("DATA")
};
}
#endif
