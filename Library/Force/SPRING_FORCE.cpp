#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Equation/MATRIX_BUNDLE.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/SPRING_FORCE.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Math/Spring_Force.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/LOG.h>
#include <Utilities/MATH.h>
#include <Utilities/OSG_HELPERS.h>
#include <iostream>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void SPRING_FORCE<TV>::
Identify_Interactions_And_Compute_Errors(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    VECTOR& right_hand_side=system.RHS(data,force,*rigid_data);
    for(int i=0;i<springs.size();i++){
        const SPRING& spring=springs[i];
        std::array<TV,2> spun_offsets,points,base_offsets,positions;
        std::array<T_SPIN,2> spins;
        for(int j=0;j<2;j++){
            auto structure=rigid_data->structures[spring.index[j]];
            spun_offsets[j]=structure->frame.orientation*spring.offset[j];
            positions[j]=structure->frame.position;
            points[j]=structure->frame*spring.offset[j];
            spins[j]=structure->twist.angular;
            base_offsets[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]).inverse()*spun_offsets[j];}
        TV relative_position=data.Minimum_Offset(points[0],points[1]);
        Spring_Force<TV>::Evaluate(spring.index,spring.stiffness,spring.target_distance,positions,spins,base_offsets,right_hand_side);
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SPRING_FORCE<TV>::
Compute_Derivatives(DATA<TV>& data,FORCE<TV>& force,MATRIX_BUNDLE<TV>& system)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    std::vector<Triplet<T>>& force_terms=system.Matrix_Block_Terms(data,force,*rigid_data);
    std::vector<Triplet<T>> terms;

    for(int i=0;i<springs.size();i++){
        const SPRING& spring=springs[i];
        std::array<TV,2> spun_offsets,points,base_offsets;
        std::array<T_SPIN,2> spins;
        for(int j=0;j<2;j++){
            auto structure=rigid_data->structures[spring.index[j]];
            spun_offsets[j]=structure->frame.orientation*spring.offset[j];
            points[j]=structure->frame*spring.offset[j];
            spins[j]=structure->twist.angular;
            base_offsets[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]).inverse()*spun_offsets[j];}
        TV relative_position=data.Minimum_Offset(points[0],points[1]);
        Spring_Force<TV>::Derivatives(spring.index,spring.stiffness,spring.target_distance,relative_position,spins,base_offsets,force_terms);
    }
    
    system.Build_Jacobian_Block(data,force,*this,*rigid_data,terms);
    system.Build_Jacobian_Block(data,force,*rigid_data,*this,terms);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SPRING_FORCE<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
    osg::Group* group=node->asGroup();
    osg::Group* relative_position_group=(osg::Group*)getNamedChild(group,Static_Name());
    if(!relative_position_group){
        relative_position_group=new osg::Group();
        relative_position_group->setName(Static_Name());
        osg::Vec4 color(0.0f,1.0f,1.0f,1.0f);
        for(int i=0;i<springs.size();i++){relative_position_group->addChild(createLine(color));}
        group->addChild(relative_position_group);}
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    for(int i=0;i<springs.size();i++){
        const SPRING& spring=springs[i];
        std::array<TV,2> points;
        for(int j=0;j<2;j++){points[j]=rigid_data->structures[spring.index[j]]->frame*spring.offset[j];}
        updateLine((osg::Geode*)relative_position_group->getChild(i),points);}
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SPRING_FORCE)
DEFINE_AND_REGISTER_PARSER(SPRING_FORCE,void)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto spring_force=simulation.force.template Find_Or_Create<SPRING_FORCE<TV>>();
    Json::Value springs=node["springs"];
    for(Json::ValueIterator it=springs.begin();it!=springs.end();it++){
        typename SPRING_FORCE<TV>::SPRING spring;
        spring.index[0]=rigid_data->Structure_Index((*it)["structure1"].asString());
        Parse_Vector((*it)["offset1"],spring.offset[0]);
        spring.index[1]=rigid_data->Structure_Index((*it)["structure2"].asString());
        Parse_Vector((*it)["offset2"],spring.offset[1]);
        spring.target_distance=(*it)["distance"].asDouble();
        spring.stiffness=(*it)["stiffness"].asDouble();
        spring_force->springs.push_back(spring);
    }
    spring_force->stored_forces.resize(0);
    return 0;
}
