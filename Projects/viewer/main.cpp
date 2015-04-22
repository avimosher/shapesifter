#include <Data/DATA.h>
#include <Data/TEST_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSE_SCENE.h>
#include <fstream>
#include <osg/Node>
#include <osg/NodeCallback>
#include <osgDB/ReadFile>
#include <osgGA/GUIEventHandler>
#include <osgGA/TrackballManipulator>
#include <osgViewer/Viewer>

using namespace Mechanics;

class AnimationHandler : public osg::NodeCallback
{
    typedef double T;
    typedef Matrix<T,3,1> TV;
public:
    int frame;
    bool animating;
    T lastTime;
    T frameTime;
    SIMULATION<TV> simulation;
    osg::Group* root;

    AnimationHandler()
        :root(new osg::Group()),animating(false),lastTime(0),frameTime((T).05)
    {
        PARSE_SCENE<TV>::Parse_Scene(std::cin,simulation);
        reset();
        root->setUpdateCallback(this);
    }

    void toggleAnimating()
    {animating=!animating;}

    void reset() {
        frame=1;
        animating=false;
        advanceFrame(0);
    }

    void advanceFrame(int increment=1) {
        if(simulation.Read(frame+increment)){frame+=increment;}
        simulation.Viewer(root);
    }

    virtual void operator()(osg::Node* node,osg::NodeVisitor* nv) {
        T time=nv->getFrameStamp()->getSimulationTime();
        if(animating && time-lastTime>frameTime){
            advanceFrame();
            lastTime=time;
        }
    }
    
};

class KeyboardEventHandler : public osgGA::GUIEventHandler
{
public:
    AnimationHandler* animation;

    KeyboardEventHandler(AnimationHandler* animation)
        :animation(animation)
    {}

    virtual bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter&) {
        switch(ea.getEventType()) {
            case(osgGA::GUIEventAdapter::KEYDOWN):{
                switch(ea.getKey()) {
                    case 's':
                        animation->advanceFrame(1);
                        return true;
                    case 19: // this is ctrl-s but I can't figure out why
                        animation->advanceFrame(-1);
                        return true;
                    case 'p':
                        animation->lastTime=ea.getTime();
                        animation->toggleAnimating();
                        return true;
                    case 'r':
                        animation->reset();
                        return true;
                    default:
                        return false;
                }}
            default:
                return false;}}

    virtual void accept(osgGA::GUIEventHandlerVisitor& v) {
        v.visit(*this);
    }
};

int main(int argc,char **argv)
{
    AnimationHandler* animation=new AnimationHandler();
    KeyboardEventHandler* keyboardEventHandler=new KeyboardEventHandler(animation);
    osgViewer::Viewer viewer;
    viewer.addEventHandler(keyboardEventHandler);
    viewer.setSceneData(animation->root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator);
    viewer.setUpViewInWindow(0,0,640,480);
    viewer.realize(); //creates windows, starts threads
    viewer.run();
    return 0;
}
