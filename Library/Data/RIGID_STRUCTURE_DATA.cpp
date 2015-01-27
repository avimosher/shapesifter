#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/OSG_HELPERS.h>
#include <osg/Geode>
#include <osg/Node>
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osgWidget/Box>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Size()
{
    return structures.size()*TWIST<TV>::STATIC_SIZE;
}
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Velocity_DOF() const
{
    return TWIST<TV>::STATIC_SIZE*structures.size();
}
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Position_DOF() const
{
    return FRAME<TV>::STATIC_SIZE*structures.size();
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Pack_Velocities(Block<Matrix<T,Dynamic,1>>& velocities)
{
    for(int i=0;i<structures.size();i++){
        velocities.template block<TWIST<TV>::STATIC_SIZE,1>(i*TWIST<TV>::STATIC_SIZE,0)=structures[i]->twist.Pack();
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities)
{
    for(int i=0;i<structures.size();i++){
        structures[i]->twist.Unpack(velocities.template block<TWIST<TV>::STATIC_SIZE,1>(i*TWIST<TV>::STATIC_SIZE,0));
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Pack_Positions(Block<Matrix<T,Dynamic,1>>& positions)
{
    for(int i=0;i<structures.size();i++){
        positions.template block<FRAME<TV>::STATIC_SIZE,1>(i*FRAME<TV>::STATIC_SIZE,0)=structures[i]->frame.Pack();
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Unpack_Positions(const Matrix<T,Dynamic,1>& positions)
{
    for(int i=0;i<structures.size();i++){
        structures[i]->frame.Unpack(positions.template block<FRAME<TV>::STATIC_SIZE,1>(i*FRAME<TV>::STATIC_SIZE,0));
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Step(const DATA<TV>& data)
{
    for(int i=0;i<structures.size();i++){
        structures[i]->frame=Updated_Frame(data,structures[i]->frame,structures[i]->twist);
        //std::cout<<i<<": "<<structures[i]->frame.position.transpose()<<std::endl;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Viewer(osg::Node* node)
{
    osg::Group* group=node->asGroup();
    osg::Group* rigid_group=NULL;
    for(int i=0;i<group->getNumChildren();i++){
        if(group->getChild(i)->getName()=="RIGID_STRUCTURE_DATA")
        {rigid_group=(osg::Group*)group->getChild(i);}}
    if(!rigid_group){
        rigid_group=new osg::Group();
        rigid_group->setName("RIGID_STRUCTURE_DATA");
        for(int i=0;i<structures.size();i++){
            auto transform=new osg::PositionAttitudeTransform();
            auto basicShapesGeode=new osg::Geode();
            if(structures[i]->collision_extent){
                auto cylinder=new osg::Cylinder(osg::Vec3(0,0,0),structures[i]->radius,2*structures[i]->collision_extent);
                auto cylinderDrawable=new osg::ShapeDrawable(cylinder);
                basicShapesGeode->addDrawable(cylinderDrawable);
                auto topSphere=new osg::Sphere(osg::Vec3(0,0,structures[i]->collision_extent),structures[i]->radius);
                auto topSphereDrawable=new osg::ShapeDrawable(topSphere);
                basicShapesGeode->addDrawable(topSphereDrawable);
                auto bottomSphere=new osg::Sphere(osg::Vec3(0,0,-structures[i]->collision_extent),structures[i]->radius);
                auto bottomSphereDrawable=new osg::ShapeDrawable(bottomSphere);
                basicShapesGeode->addDrawable(bottomSphereDrawable);
            }
            else{
                auto unitSphere=new osg::Sphere(osg::Vec3(0,0,0),structures[i]->radius);
                auto unitSphereDrawable=new osg::ShapeDrawable(unitSphere);
                basicShapesGeode->addDrawable(unitSphereDrawable);
            }
            transform->addChild(basicShapesGeode);
            rigid_group->addChild(transform);
        }
        group->addChild(rigid_group);
    }
    for(int i=0;i<structures.size();i++){
        auto transform=(osg::PositionAttitudeTransform*)rigid_group->getChild(i);
        OSG_HELPERS<TV>::Initialize_Transform(structures[i]->frame,transform);
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RIGID_STRUCTURE_DATA)
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_DATA)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE_DATA,void)
{
    simulation.data.insert({RIGID_STRUCTURE_DATA<TV>::Static_Name(),std::make_shared<RIGID_STRUCTURE_DATA<TV>>()});
    return 0;
}
