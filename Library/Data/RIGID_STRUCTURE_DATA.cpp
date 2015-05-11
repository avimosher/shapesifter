#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/OSG_HELPERS.h>
#include <osg/Geode>
#include <osg/Node>
#include <osg/PolygonMode>
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osgWidget/Box>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Size()
{
    return structures.size();
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
        TV p1=structures[i]->frame.position;
        structures[i]->frame=Updated_Frame(data,structures[i]->frame,structures[i]->twist);
        //std::cout<<"Kick in center of mass for "<<i<<": "<<(structures[i]->frame.position-p1).transpose()<<std::endl;
        //std::cout<<i<<": "<<structures[i]->frame.position.transpose()<<std::endl;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Inertia(const T dt,std::vector<Triplet<T>>& force_terms,Matrix<T,Dynamic,1>& rhs)
{
    T one_over_dt=1/dt;
    for(int i=0;i<structures.size();i++){
        DiagonalMatrix<T,s> inertia_matrix=one_over_dt*structures[i]->Inertia_Matrix();
        Flatten_Term(i,i,inertia_matrix,force_terms);
        //std::cout<<"Twist "<<i<<": "<<structures[i]->twist.Pack().transpose()<<std::endl;
        rhs.template block<s,1>(s*i,0)=-(inertia_matrix*structures[i]->twist.Pack());
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
                cylinderDrawable->setColor(osg::Vec4(1.0f,1.0f,1.0f,0.5f));
                basicShapesGeode->addDrawable(cylinderDrawable);
                auto topSphere=new osg::Sphere(osg::Vec3(0,0,structures[i]->collision_extent),structures[i]->radius);
                auto topSphereDrawable=new osg::ShapeDrawable(topSphere);
                topSphereDrawable->setColor(osg::Vec4(1.0f,1.0f,1.0f,0.5f));
                basicShapesGeode->addDrawable(topSphereDrawable);
                auto bottomSphere=new osg::Sphere(osg::Vec3(0,0,-structures[i]->collision_extent),structures[i]->radius);
                auto bottomSphereDrawable=new osg::ShapeDrawable(bottomSphere);
                bottomSphereDrawable->setColor(osg::Vec4(1.0f,1.0f,1.0f,0.5f));
                basicShapesGeode->addDrawable(bottomSphereDrawable);
            }
            else{
                auto unitSphere=new osg::Sphere(osg::Vec3(0,0,0),structures[i]->radius);
                auto unitSphereDrawable=new osg::ShapeDrawable(unitSphere);
                unitSphereDrawable->setColor(osg::Vec4(1.0f,1.0f,1.0f,0.5f));
                basicShapesGeode->addDrawable(unitSphereDrawable);
            }
            /*auto program=new osg::Program;
            auto fragmentShader=new osg::Shader(osg::Shader::FRAGMENT);
            fragmentShader->setShaderSource(
                "in vec2 vTexCoord;\n"
                "void main(void)\n"
                "{\n"
                "    float x;x=fract(vTexCoord.x);\n"
                "    if(x>0.5){gl_FragColor=vec4(1.0,1.0,1.0,1.0);}\n"
                "    else{gl_FragColor=vec4(0.0,0.0,0.0,0.0);}\n"
//                "    gl_FragColor=gl_FragCoord;\n"
//                "    gl_FragColor=color*texture2D(base,gl_TexCoord[0]);\n"
//                "    gl_FragColor=ColorFilter(gl_FragColor);\n"
                "}\n"
            );
            fragmentShader->setName("fragment");
            program->addShader(fragmentShader);
            auto vertexShader=new osg::Shader(osg::Shader::VERTEX);
            vertexShader->setShaderSource(
                "in vec2 aTexCoord;\n"
                "out vec2 vTexCoord;\n"
                "void main(void)\n"
                "{\n"
                "    vTexCoord=aTexCoord;\n"
                "}\n"
            );
            vertexShader->setName("vertex");
            program->addShader(vertexShader);
            auto stateSet=basicShapesGeode->getOrCreateStateSet();
            stateSet->setAttributeAndModes(program,osg::StateAttribute::ON);*/
            basicShapesGeode->getOrCreateStateSet()->setAttribute(new osg::PolygonMode(osg::PolygonMode::FRONT_AND_BACK,osg::PolygonMode::LINE));
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
    simulation.data.push_back(std::make_shared<RIGID_STRUCTURE_DATA<TV>>());
    return 0;
}
///////////////////////////////////////////////////////////////////////
namespace Eigen{
namespace internal{
template<> AlignedBox<double,3> bounding_box(const std::shared_ptr<RIGID_STRUCTURE<Matrix<double,3,1>>> structure)
{
    typedef double T;
    const static int d=3;
    Matrix<T,d,1> minimum=structure->frame*(-Matrix<T,d,1>::Unit(2)*structure->collision_extent);
    Matrix<T,d,1> maximum=structure->frame*(Matrix<T,d,1>::Unit(2)*structure->collision_extent);
    return AlignedBox<T,3>(minimum.cwiseMin(maximum)-Matrix<double,3,1>::Constant(structure->radius),minimum.cwiseMax(maximum)+Matrix<double,3,1>::Constant(structure->radius));
}
template<> AlignedBox<float,3> bounding_box(const std::shared_ptr<RIGID_STRUCTURE<Matrix<float,3,1>>> structure)
{
    typedef float T;
    const static int d=3;
    Matrix<T,d,1> minimum=structure->frame*(-Matrix<T,d,1>::Unit(2)*structure->collision_extent);
    Matrix<T,d,1> maximum=structure->frame*(Matrix<T,d,1>::Unit(2)*structure->collision_extent);
    return AlignedBox<T,3>(minimum.cwiseMin(maximum)-Matrix<float,3,1>::Constant(structure->radius),minimum.cwiseMax(maximum)+Matrix<float,3,1>::Constant(structure->radius));
}
}
}
