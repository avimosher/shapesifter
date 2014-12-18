#include <Data/DATA.h>
#include <Data/TEST_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSE_SCENE.h>
#include <fstream>
#include <osg/Node>
#include <osgDB/ReadFile>
#include <osgGA/GUIEventHandler>
#include <osgGA/TrackballManipulator>
#include <osgViewer/Viewer>

using namespace Mechanics;

class KeyboardEventHandler : public osgGA::GUIEventHandler
{
    typedef double T;
    typedef Matrix<T,3,1> TV;
    int frame;
public:
    SIMULATION<TV> simulation;
    osg::Group* root;

    KeyboardEventHandler()
        :frame(1),root(new osg::Group()) 
    {
        std::ifstream test_config("brownian.json",std::ifstream::in);
        PARSE_SCENE<TV>::Parse_Scene(test_config,simulation);
        simulation.Read(frame);
        simulation.data.Viewer(root);
    }

    virtual bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter&) {
        switch(ea.getEventType()) {
            case(osgGA::GUIEventAdapter::KEYDOWN):{
                switch(ea.getKey()) {
                    case 's':
                        simulation.Read(++frame);
                        simulation.data.Viewer(root);
                        return true;
                        break;
                    default:
                        return false;
                }}
            default:
                return false;}}

    virtual void accept(osgGA::GUIEventHandlerVisitor& v) {
        v.visit(*this);
    }
};

int main()
{
    KeyboardEventHandler* keyboardEventHandler=new KeyboardEventHandler();
    osgViewer::Viewer viewer;
    viewer.addEventHandler(keyboardEventHandler);
    viewer.setSceneData(keyboardEventHandler->root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator);
    viewer.realize(); //creates windows, starts threads
    while(!viewer.done()){viewer.frame();}

    return 0;
}
