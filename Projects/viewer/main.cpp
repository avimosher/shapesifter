#include <Data/DATA.h>
#include <Data/TEST_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSE_SCENE.h>
#include <fstream>
#include <iomanip>
#include <osg/Node>
#include <osg/NodeCallback>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgGA/GUIEventHandler>
#include <osgGA/TrackballManipulator>
#include <osgViewer/Renderer>
#include <osgViewer/Viewer>

using namespace Mechanics;

class WindowCaptureCallback:public osg::Camera::DrawCallback
{
  public:
    WindowCaptureCallback(GLenum readBuffer)
      :_readBuffer(readBuffer)
    {
      _image=new osg::Image;
    }

    virtual void operator()(osg::RenderInfo& renderInfo) const
    {
      glReadBuffer(_readBuffer);
      OpenThreads::ScopedLock<OpenThreads::Mutex> lock(_mutex);
      osg::GraphicsContext* gc=renderInfo.getState()->getGraphicsContext();
      if(gc->getTraits()){
        GLenum pixelFormat;
        if(gc->getTraits()->alpha){pixelFormat=GL_RGBA;}
        else{pixelFormat=GL_RGB;}
        
        int width=gc->getTraits()->width;
        int height=gc->getTraits()->height;
        _image->readPixels(0,0,width,height,pixelFormat,GL_UNSIGNED_BYTE);
      }

      if(!_fileName.empty()){
        osgDB::writeImageFile(*_image,_fileName);
      }
    }

    void setFrame(const std::string& name)
    {_fileName=name;}

  protected:
    GLenum _readBuffer;
    std::string _fileName;
    osg::ref_ptr<osg::Image> _image;
    mutable OpenThreads::Mutex _mutex;
};

class AnimationHandler : public osg::NodeCallback
{
    typedef double T;
    typedef Matrix<T,3,1> TV;
public:
    int frame;
    bool animating;
    bool writing;
    T lastTime;
    T frameTime;
    SIMULATION<TV> simulation;
    osg::Group* root;
    WindowCaptureCallback *capture;
    osgViewer::Viewer* pViewer;

    AnimationHandler(osgViewer::Viewer* viewer)
        :root(new osg::Group()),animating(false),writing(false),lastTime(0),frameTime((T).05),pViewer(viewer)
    {
        PARSE_SCENE<TV>::Parse_Scene(std::cin,simulation);
        reset();
        root->setUpdateCallback(this);
        capture=new WindowCaptureCallback(GL_FRONT);
    }

    void toggleAnimating()
    {animating=!animating;}

    void toggleWriting()
    {writing=!writing;
      if(writing){pViewer->getCamera()->setPostDrawCallback(capture);}
      else{pViewer->getCamera()->setPostDrawCallback(NULL);}
    }

    void reset() {
        frame=1;
        animating=false;
        writing=false;
        advanceFrame(0);
    }

    void advanceFrame(int increment=1) {
        if(simulation.Read(frame+increment)){frame+=increment;}
        if(writing){
          std::stringstream stream;
          stream<<simulation.output_directory+"/image."<<std::setw(5)<<std::setfill('0')<<frame<<".png";
          capture->setFrame(stream.str());}
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
                    case 'w':
                        animation->toggleWriting();
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
    osgViewer::Viewer viewer;
    AnimationHandler* animation=new AnimationHandler(&viewer);
    KeyboardEventHandler* keyboardEventHandler=new KeyboardEventHandler(animation);
    viewer.addEventHandler(keyboardEventHandler);
    viewer.setSceneData(animation->root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator);
    viewer.setUpViewInWindow(0,0,640,480);
    viewer.realize(); //creates windows, starts threads
    viewer.run();
    return 0;
}
