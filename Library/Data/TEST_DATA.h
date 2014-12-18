//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class TEST_DATA
///////////////////////////////////////////////////////////////////////
#ifndef __TEST_DATA__
#define __TEST_DATA__

#include <Data/DATA_TYPE.h>
#include <iostream>
#include <cereal/archives/binary.hpp>
#include <osg/Geode>
#include <osg/Node>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>
#include <osgWidget/Box>

namespace Mechanics{

template<class TV>
class TEST_DATA:public DATA_TYPE<TV>
{
    typedef typename TV::Scalar T;

public:
    T internal_data; // v. simple

    TEST_DATA(){internal_data=5;}
    ~TEST_DATA(){}

    virtual int Size(){
        return 1;
    }

    virtual Matrix<T,Dynamic,1> Variables() {
        Matrix<T,Dynamic,1> variables(1,1);
        variables(0,0)=internal_data;
        return variables;
    }

    virtual void Step(const Matrix<T,Dynamic,1>& variables) {
        internal_data=variables(0,0);
    }

    template<class Archive>
    void serialize(Archive& archive) {
        std::cout<<"Serialize"<<std::endl;
        archive(cereal::make_nvp("internal_data",internal_data));
    }
    
    virtual T Print() {
        return internal_data;
    }

    virtual void Viewer(osg::Node* node) {
        std::cout<<"Viewing"<<std::endl;
        osg::Group* group=node->asGroup();
        osg::PositionAttitudeTransform* transform=NULL;
        for(int i=0;i<group->getNumChildren();i++){
            if(group->getChild(i)->getName()=="TEST_DATA") {
                transform=(osg::PositionAttitudeTransform*)group->getChild(i);
                break;
            }
        }
        std::cout<<transform<<std::endl;
        if(!transform){
            transform=new osg::PositionAttitudeTransform();
            osg::Geode* basicShapesGeode=new osg::Geode();
            osg::Box* unitCube=new osg::Box(osg::Vec3(0,0,0),1.0f);
            osg::ShapeDrawable* unitCubeDrawable=new osg::ShapeDrawable(unitCube);
            basicShapesGeode->addDrawable(unitCubeDrawable);
            transform->addChild(basicShapesGeode);
            transform->setName("TEST_DATA");
            group->addChild(transform);
        }
        osg::Vec3 pos(internal_data,internal_data,internal_data);
        transform->setPosition(pos);
    }

    DEFINE_TYPE_NAME("TEST_DATA")
};
}
#endif
