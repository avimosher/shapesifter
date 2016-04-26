#ifndef __OSG_HELPERS__
#define __OSG_HELPERS__

#include <Data/FRAME.h>
#include <Eigen/Geometry>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/PositionAttitudeTransform>
#include <functional>
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

inline osg::Vec4 colorMap(const std::string& name)
{
    static const std::vector<osg::Vec4> colors={
        osg::Vec4(0.271,0.647,0.424,1),
        osg::Vec4(0.627,0.863,0.722,1),
        osg::Vec4(0.424,0.753,0.557,1),
        osg::Vec4(0.165,0.561,0.325,1),
        osg::Vec4(0.078,0.455,0.231,1),
        osg::Vec4(0.251,0.471,0.569,1),
        osg::Vec4(0.608,0.753,0.82,1),
        osg::Vec4(0.396,0.588,0.675,1),
        osg::Vec4(0.161,0.392,0.49,1),
        osg::Vec4(0.086,0.306,0.4,1),
        osg::Vec4(0.902,0.655,0.376,1),
        osg::Vec4(1,0.871,0.729,1),
        osg::Vec4(1,0.796,0.565,1),
        osg::Vec4(0.776,0.518,0.227,1),
        osg::Vec4(0.635,0.388,0.11,1),
        osg::Vec4(0.902,0.478,0.376,1),
        osg::Vec4(1,0.78,0.729,1),
        osg::Vec4(1,0.647,0.565,1),
        osg::Vec4(0.776,0.333,0.227,1)};
//        osg::Vec4(0.635,0.212,0.11,1)};
    return colors[std::hash<std::string>()(name)%colors.size()];
}

inline osg::Geode* createLine(osg::Vec4 color)
{
    auto errorLineGeometry=new osg::Geometry();
    auto errorVertices=new osg::Vec3Array(2);
    errorLineGeometry->setVertexArray(errorVertices);
    auto errorColors=new osg::Vec4Array;
    errorColors->push_back(color);
    errorLineGeometry->setColorArray(errorColors);
    errorLineGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
    auto errorNormals=new osg::Vec3Array;
    errorNormals->push_back(osg::Vec3f(0.0f,-1.0f,0.0f));
    errorLineGeometry->setNormalArray(errorNormals);
    errorLineGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);

    errorLineGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,2));

    osg::StateSet* errorStateset=new osg::StateSet;
    errorStateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
    errorLineGeometry->setStateSet(errorStateset);

    auto errorLineGeode=new osg::Geode();
    errorLineGeode->addDrawable(errorLineGeometry);
    return errorLineGeode;
}

template<class TV>
inline void updateLine(osg::Geode* lineGeode,const std::array<TV,2>& points)
{
    auto lineGeometry=(osg::Geometry*)lineGeode->getDrawable(0);
    auto vertices=(osg::Vec3Array*)lineGeometry->getVertexArray();
    for(int i=0;i<points.size();i++){
        (*vertices)[i].set(points[i](0),points[i](1),points[i](2));}
    lineGeometry->setVertexArray(vertices);
}

}

#endif
