#include <Data/DATA.h>
#include <Data/TEST_DATA.h>
#include <osg/Node>
#include <osgWidget/Box>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#include <osgGA/GUIEventHandler>
#include <osg/PositionAttitudeTransform>
#include <osgDB/ReadFile>
#include <osgGA/TrackballManipulator>
#include <osgViewer/Viewer>

using namespace Mechanics;
osg::PositionAttitudeTransform* transform;

class KeyboardEventHandler : public osgGA::GUIEventHandler
{
    typedef double T;
    typedef Matrix<T,1,1> TV;
    DATA<TV> data;
    TEST_DATA<TV>* test_data;
    int frame;
public:
    KeyboardEventHandler() {
        frame=1;
        test_data=new TEST_DATA<TV>();
        data.push_back(std::unique_ptr<DATA_TYPE<TV>>(test_data));
        data.Read(frame);
        osg::Vec3 pos(data.Print_All(),0,0);
        transform->setPosition(pos);
    }

    virtual bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter&) {
        switch(ea.getEventType()) {
            case(osgGA::GUIEventAdapter::KEYDOWN):
                {
                    switch(ea.getKey()) {
                        case 's':
                            {
                            data.Read(++frame);
                            osg::Vec3 pos(data.Print_All(),0,0);
                            transform->setPosition(pos);
                            std::cout<<"s pressed"<<std::endl;
                            std::cout<<"Stuff"<<std::endl;
                            return true;
                            }
                            break;
                        default:
                            return false;
                    }
                }
            default:
                return false;
        }
                        
    }

    virtual void accept(osgGA::GUIEventHandlerVisitor& v) {
        v.visit(*this);
    }
};

int main()
{
    

    osg::Node* cessnaNode=NULL;
    osg::Node* tankNode=NULL;
    osg::Box* unitCube=new osg::Box(osg::Vec3(0,0,0),1.0f);

    osg::ShapeDrawable* unitCubeDrawable=new osg::ShapeDrawable(unitCube);
    osg::Geode* basicShapesGeode=new osg::Geode();
    basicShapesGeode->addDrawable(unitCubeDrawable);
    transform=new osg::PositionAttitudeTransform();

    osg::Group* root=new osg::Group();
    transform->addChild(basicShapesGeode);
    root->addChild(transform);
    
    /*osg::PositionAttitudeTransform* tankXform=new osg::PositionAttitudeTransform();
    root->addChild(tankXform);
    tankXform->addChild(tankNode);*/
    
    KeyboardEventHandler* keyboardEventHandler=new KeyboardEventHandler();
    osgViewer::Viewer viewer;
    viewer.addEventHandler(keyboardEventHandler);
    viewer.setSceneData(root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator);
    viewer.realize(); //creates windows, starts threads
    
    int count=0;

    while(!viewer.done()){
        count++;
        osg::Vec3 tankPosit(count/(float)40,0,0);
        //tankXform->setPosition(tankPosit);
        viewer.frame();
    }

    return 0;
}
