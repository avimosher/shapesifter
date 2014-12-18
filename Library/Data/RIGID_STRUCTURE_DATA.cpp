//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE_DATA
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <osg/Geode>
#include <osg/Node>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>
#include <osgWidget/Box>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> RIGID_STRUCTURE_DATA<TV>::
RIGID_STRUCTURE_DATA()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Size()
{
    return structures.size()*RIGID_STRUCTURE<TV>::STATIC_SIZE;
}
///////////////////////////////////////////////////////////////////////
template<class TV> Matrix<typename TV::Scalar,Dynamic,1> RIGID_STRUCTURE_DATA<TV>::
Variables()
{
    Matrix<T,Dynamic,1> packed(RIGID_STRUCTURE<TV>::STATIC_SIZE*structures.size(),1);
    for(int i=0;i<structures.size();i++){
        packed.template block<RIGID_STRUCTURE<TV>::STATIC_SIZE,1>(i*RIGID_STRUCTURE<TV>::STATIC_SIZE,0)=structures[i]->frame.Pack();
    }
    return packed;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Step(const Matrix<T,Dynamic,1>& variables)
{
    for(int i=0;i<structures.size();i++){
        structures[i]->frame.Unpack(variables.template block<RIGID_STRUCTURE<TV>::STATIC_SIZE,1>(i*RIGID_STRUCTURE<TV>::STATIC_SIZE,0));
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Viewer(osg::Node* node)
{
    osg::Group* group=node->asGroup();
    osg::Group* rigid_group=NULL;
    for(int i=0;i<group->getNumChildren();i++){
        if(group->getChild(i)->getName()=="RIGID_STRUCTURE_DATA"){
            rigid_group=(osg::Group*)group->getChild(i);
        }
    }
    if(!rigid_group){
        rigid_group=new osg::Group();
        rigid_group->setName("RIGID_STRUCTURE_DATA");
        for(int i=0;i<structures.size();i++){
            auto transform=new osg::PositionAttitudeTransform();
            auto unitCube=new osg::Box(osg::Vec3(0,0,0),1.0f);
            auto unitCubeDrawable=new osg::ShapeDrawable(unitCube);
            auto basicShapesGeode=new osg::Geode();
            basicShapesGeode->addDrawable(unitCubeDrawable);
            transform->addChild(basicShapesGeode);
            rigid_group->addChild(transform);
        }
        group->addChild(rigid_group);
    }
    for(int i=0;i<structures.size();i++){
        auto transform=(osg::PositionAttitudeTransform*)rigid_group->getChild(i);
        osg::Vec3 pos;
        for(int j=0;j<3;j++){pos[j]=structures[i]->frame.position[j];} // TODO: won't work for 1d, 2d
        transform->setPosition(pos);
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RIGID_STRUCTURE_DATA)
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_DATA)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE_DATA)
{
    simulation.data.insert({RIGID_STRUCTURE_DATA<TV>::Static_Name(),std::make_shared<RIGID_STRUCTURE_DATA<TV>>()});
}
