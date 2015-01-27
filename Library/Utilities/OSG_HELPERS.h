#ifndef __OSG_HELPERS__
#define __OSG_HELPERS__

#include <Data/FRAME.h>
#include <Eigen/Geometry>
#include <osg/Geode>
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
        osg::Quat quat(frame.orientation.x(),frame.orientation.y(),frame.orientation.z(),frame.orientation.w());
        transform->setAttitude(quat);
    }
};

inline osg::Node* getNamedChild(osg::Group* group,const std::string& name)
{
    for(int i=0;i<group->getNumChildren();i++){
        if(group->getChild(i)->getName()==name){
            return group->getChild(i);
        }}
    return 0;
}

inline osg::Drawable* getNamedDrawable(osg::Geode* geode,const std::string& name)
{
    for(int i=0;i<geode->getNumDrawables();i++){
        if(geode->getDrawable(i)->getName()==name){
            return geode->getDrawable(i);
        }}
    return 0;
}

}

#endif
