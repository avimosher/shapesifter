#include <osg/Node>
#include <osg/PositionAttitudeTransform>
#include <osgDB/ReadFile>
#include <osgGA/TrackballManipulator>
#include <osgViewer/Viewer>


int main()
{
    osg::Node* cessnaNode=NULL;
    osg::Node* tankNode=NULL;
    cessnaNode=osgDB::readNodeFile("cessna.osg");
    tankNode=osgDB::readNodeFile("spaceship.osgt");

    osg::Group* root=new osg::Group();
    root->addChild(cessnaNode);
    
    osg::PositionAttitudeTransform* tankXform=new osg::PositionAttitudeTransform();
    root->addChild(tankXform);
    tankXform->addChild(tankNode);
    
    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator);
    viewer.realize(); //creates windows, starts threads
    
    int count=0;

    while(!viewer.done()){
        count++;
        osg::Vec3 tankPosit(count/(float)40,0,0);
        tankXform->setPosition(tankPosit);
        viewer.frame();
    }

    return 0;
}
