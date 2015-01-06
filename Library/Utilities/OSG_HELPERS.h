#ifndef __OSG_HELPERS__
#define __OSG_HELPERS__

#include <Data/FRAME.h>
#include <Eigen/Geometry>
#include <osg/PositionAttitudeTransform>
#include <iostream>

namespace Mechanics{

template<class TV>
class OSG_HELPERS
{
public:
    static void Initialize_Transform(const FRAME<TV>& frame,osg::PositionAttitudeTransform* transform)
    {}
};


template<class T>
class OSG_HELPERS<Matrix<T,3,1>>
{
  public:
    static void Initialize_Transform(const FRAME<Matrix<T,3,1>>& frame,osg::PositionAttitudeTransform* transform)
    {
        osg::Vec3 pos;
        for(int j=0;j<3;j++){pos[j]=frame.position[j];}
        transform->setPosition(pos);
//      osg::Vec3 axis;axis[0]=1;axis[1]=1;axis[2]=1;
//      osg::Quat quat(M_PI/4,axis);
        //AngleAxis angleAxis(frame.orientation);
//std::cout<<"Quaternion: "<<std::endl<<frame.orientation.coeffs()<<std::endl;
//transform->setAttitude(quat);
osg::Quat quat(frame.orientation.x(),frame.orientation.y(),frame.orientation.z(),frame.orientation.w());
transform->setAttitude(quat);
    }
};

}

#endif
