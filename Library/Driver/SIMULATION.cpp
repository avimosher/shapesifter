#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Evolution/EVOLUTION.h>
#include <Force/FORCE.h>
#include <Force/FORCE_TYPE.h>
#include <Utilities/OSG_HELPERS.h>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <Eigen/Geometry>
#include <osg/Geode>
#include <osg/MatrixTransform>
#include <osg/Projection>
#include <osg/ShapeDrawable>
#include <osgText/Text>
#include <sys/stat.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> SIMULATION<TV>::
SIMULATION()
    :data(*new DATA<TV>()),evolution(*new EVOLUTION<TV>()),force(*new FORCE<TV>()),
    current_frame(0),restart_frame(0),output_number(0),time(0),restart(false),substeps(false),output_directory("."),title("")
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> SIMULATION<TV>::
~SIMULATION()
{
    delete &data;
    delete &evolution;
    delete &force;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SIMULATION<TV>::
Write(const std::string& frame_title)
{
    title=frame_title;
    std::ostringstream stringStream;
    mkdir(output_directory.c_str(),0777);
    stringStream<<output_directory<<"/frame."<<output_number++;
    std::ofstream output(stringStream.str().c_str(),std::ios::out);
    cereal::BinaryOutputArchive archive(output);
    archive(time,title,data,force);
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> bool SIMULATION<TV>::
Read(const int frame)
{
    std::ostringstream stringStream;
    stringStream<<output_directory<<"/frame."<<frame;
    std::ifstream input(stringStream.str().c_str(),std::ios::in);
    if(!input.is_open()){return false;}
    cereal::BinaryInputArchive archive(input);
    archive(time,title,data,force);
    return true;
}
/////////////////////////////////////////////////////////////////////// 
template<class TV> void SIMULATION<TV>::
Viewer(osg::Group*& root)
{
    root=root?root:new osg::Group();

    auto hud_projection_matrix=(osg::Projection*)getNamedChild(root,"HUD");
    if(!hud_projection_matrix){
        hud_projection_matrix=new osg::Projection;
        hud_projection_matrix->setName("HUD");
        hud_projection_matrix->setMatrix(osg::Matrix::ortho2D(0,1024,0,768));
        auto hud_modelview_matrix=new osg::MatrixTransform;
        hud_modelview_matrix->setMatrix(osg::Matrix::identity());
        hud_modelview_matrix->setReferenceFrame(osg::Transform::ABSOLUTE_RF);
        root->addChild(hud_projection_matrix);
        hud_projection_matrix->addChild(hud_modelview_matrix);
        auto hud_geode=new osg::Geode();
        hud_modelview_matrix->addChild(hud_geode);

        auto hud_background_geometry=new osg::Geometry();
        auto hud_background_vertices=new osg::Vec3Array();
        hud_background_vertices->push_back(osg::Vec3(0,0,-1));
        hud_background_vertices->push_back(osg::Vec3(1024,0,-1));
        hud_background_vertices->push_back(osg::Vec3(1024,200,-1));
        hud_background_vertices->push_back(osg::Vec3(0,200,-1));

        auto hud_background_indices=new osg::DrawElementsUInt(osg::PrimitiveSet::POLYGON,0);
        for(int i=0;i<4;i++){hud_background_indices->push_back(i);}
        auto hud_colors=new osg::Vec4Array();
        hud_colors->push_back(osg::Vec4(0.8f,0.8f,0.8f,0.8f));

        auto texcoords=new osg::Vec2Array(4);
        (*texcoords)[0].set(0.0f,0.0f);
        (*texcoords)[1].set(1.0f,0.0f);
        (*texcoords)[2].set(1.0f,1.0f);
        (*texcoords)[3].set(0.0f,1.0f);
        hud_background_geometry->setTexCoordArray(0,texcoords);
        auto hud_normals=new osg::Vec3Array;
        hud_normals->push_back(osg::Vec3(0.0f,0.0f,1.0f));
        hud_background_geometry->setNormalArray(hud_normals);
        hud_background_geometry->setNormalBinding(osg::Geometry::BIND_OVERALL);
        hud_background_geometry->addPrimitiveSet(hud_background_indices);
        hud_background_geometry->setVertexArray(hud_background_vertices);
        hud_background_geometry->setColorArray(hud_colors);
        hud_background_geometry->setColorBinding(osg::Geometry::BIND_OVERALL);
        hud_geode->addDrawable(hud_background_geometry);

        auto hud_state_set=new osg::StateSet();
        hud_geode->setStateSet(hud_state_set);
        hud_state_set->setMode(GL_BLEND,osg::StateAttribute::ON);
        hud_state_set->setMode(GL_DEPTH_TEST,osg::StateAttribute::OFF);
        hud_state_set->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);

        auto text_one=new osgText::Text();
        hud_geode->addDrawable(text_one);
        text_one->setName("title");
        text_one->setCharacterSize(25);
        text_one->setText(title);
        text_one->setAxisAlignment(osgText::Text::SCREEN);
        text_one->setPosition(osg::Vec3(360,165,-1.5));
        text_one->setColor(osg::Vec4(199,77,15,1));}

    auto hud_modelview_matrix=(osg::MatrixTransform*)hud_projection_matrix->getChild(0);
    auto hud_geode=(osg::Geode*)hud_modelview_matrix->getChild(0);
    auto text_one=(osgText::Text*)getNamedDrawable(hud_geode,"title");
    text_one->setText(title);
    
    data.Viewer(root);
    force.Viewer(data,root);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SIMULATION)
