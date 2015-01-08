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

    AnimationHandler(const std::string &scenefile)
        :frame(1),root(new osg::Group()),animating(false),lastTime(0),frameTime((T).1)
    {
        std::ifstream test_config(scenefile,std::ifstream::in);
        PARSE_SCENE<TV>::Parse_Scene(test_config,simulation);
        simulation.Read(frame);
        simulation.data.Viewer(root);
        root->setUpdateCallback(this);
    }

    void toggleAnimating()
    {animating=!animating;}

    void advanceFrame(int increment=1) {
        if(simulation.Read(frame+increment)){frame+=increment;}
        simulation.data.Viewer(root);
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
                {
                    //int increment=1;
                    //std::cout<<ea.getModKeyMask()<<" "<<osgGA::GUIEventAdapter::MODKEY_CTRL<<std::endl;
                    //if(ea.getModKeyMask() & osgGA::GUIEventAdapter::MODKEY_CTRL){increment=-1;}
                    animation->advanceFrame(1);
                    return true;
                }
                case 19: // this is ctrl-s but I can't figure out why
                    animation->advanceFrame(-1);
                    return true;
                case 'p':
                    animation->lastTime=ea.getTime();
                    animation->toggleAnimating();
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
    AnimationHandler* animation=new AnimationHandler(argv[1]);
    KeyboardEventHandler* keyboardEventHandler=new KeyboardEventHandler(animation);
    osgViewer::Viewer viewer;
    viewer.addEventHandler(keyboardEventHandler);
    viewer.setSceneData(animation->root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator);
    viewer.realize(); //creates windows, starts threads
    viewer.run();
    //while(!viewer.done()){viewer.frame();}

    return 0;
}
