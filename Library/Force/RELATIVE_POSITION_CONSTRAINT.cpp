#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/RELATIVE_POSITION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/MATH.h>
#include <Utilities/OSG_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <osg/Geometry>
#include <math.h>
#include <osg/Geode>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.Find("RIGID_STRUCTURE_DATA"));
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
    constraint_rhs.resize(constraints.size(),1);
    for(int i=0;i<constraints.size();i++){
        const CONSTRAINT& constraint=constraints[i];
        int body_index1=constraint.s1;
        int body_index2=constraint.s2;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        FRAME<TV> frame1=rigid_structure1->frame;
        FRAME<TV> frame2=rigid_structure2->frame;
        TV direction=data.Minimum_Offset(frame1*constraint.v1,frame2*constraint.v2);
        T distance=direction.norm();
        std::cout<<"Constraint distance: "<<distance<<std::endl;
        Safe_Normalize(direction);
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index2,direction.transpose()*index_map.Velocity_Map(*rigid_structure2,constraint.v2)));
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index1,-direction.transpose()*index_map.Velocity_Map(*rigid_structure1,constraint.v1)));
        constraint_rhs(i,0)=constraint.target_distance-distance;
    }
    constraint_terms.resize(constraints.size(),RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE*rigid_data->structures.size());
    Flatten_Matrix(terms,constraint_terms);
    stored_forces.resize(constraints.size());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Special(DATA<TV>& data,const T dt,const T target_time)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.Find("RIGID_STRUCTURE_DATA"));
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
    for(int i=0;i<constraints.size();i++){
        const CONSTRAINT& constraint=constraints[i];
        int body_index1=constraint.s1;
        int body_index2=constraint.s2;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        FRAME<TV> frame1=rigid_structure1->frame;
        FRAME<TV> frame2=rigid_structure2->frame;
        TV x1=frame1*constraint.v1;
        TV x2=frame2*constraint.v2;
        TV direction=data.Minimum_Offset(x1,x2);
        TV offset1=frame1*constraint.v1; // this is not actually right.  I need the old values for rotation.
        TV offset2=frame2*constraint.v2;

        auto a=rigid_structure1->twist.angular;
        auto norm_a=a.norm();
        auto half_data_da=a; // zero for translation, componentwise values for rotation
        auto dw_da=-sin(norm_a/2)/(norm_a)*half_data_da.transpose();
        auto dq_da=cos(norm_a/2)*a*half_data_da.transpose()/(sqr(norm_a))+sin(norm_a/2)*(Matrix<T,t,t>::Identity()/norm_a-a*half_data_da.transpose()/cube(norm_a));
        auto w=cos(norm_a/2);
        auto q=sin(norm_a/2)*a/norm_a;
        Matrix<T,d,t+d> dvt_da; // identity for translation parts, zero for rotation
        dvt_da.template block<d,d>(0,0).setIdentity();
        dvt_da.template block<d,t>(0,d).setZero();

        auto dvror_da_partial=-2*Cross_Product_Matrix(offset2)*(q*dw_da+w*dq_da);
        auto dx2_da=dvt_da;//+dvror_da;
        dx2_da.template block<d,t>(0,d)=dvror_da_partial;
        auto dx1_da=dx2_da;
        //auto dx1_da=dvt_da+dvror_da; // TODO: specialize these to the right unknowns
        auto dx2mx1_da=dx2_da; // TODO: which a?
        auto dd_da=(dx2_da-dx1_da)/norm_a-direction*(2*(x2-x1).transpose()*dx2mx1_da)/cube(norm_a);
        auto final=direction.transpose()*dd_da+direction.transpose()*(dx2_da-dx1_da);
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
    osg::Group* group=node->asGroup();
    osg::Group* relative_position_group=(osg::Group*)getNamedChild(group,Static_Name());
    if(!relative_position_group){
        relative_position_group=new osg::Group();
        relative_position_group->setName(Static_Name());
        for(int i=0;i<constraints.size();i++){
            auto lineGeometry=new osg::Geometry();
            auto vertices=new osg::Vec3Array(2);
            lineGeometry->setVertexArray(vertices);
            auto colors=new osg::Vec4Array;
            colors->push_back(osg::Vec4(1.0f,1.0f,0.0f,1.0f));
            lineGeometry->setColorArray(colors);
            lineGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
            auto normals=new osg::Vec3Array;
            normals->push_back(osg::Vec3f(0.0f,-1.0f,0.0f));
            lineGeometry->setNormalArray(normals);
            lineGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);

            lineGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,2));
            auto lineGeode=new osg::Geode();
            lineGeode->addDrawable(lineGeometry);
            relative_position_group->addChild(lineGeode);
        }
        group->addChild(relative_position_group);
    }
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.Find("RIGID_STRUCTURE_DATA"));
    for(int i=0;i<constraints.size();i++){
        auto lineGeode=(osg::Geode*)relative_position_group->getChild(i);
        auto lineGeometry=(osg::Geometry*)lineGeode->getDrawable(0);
        auto vertices=(osg::Vec3Array*)lineGeometry->getVertexArray();
        const CONSTRAINT& constraint=constraints[i];
        int body_index1=constraint.s1;
        int body_index2=constraint.s2;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        auto firstAttachment=rigid_structure1->frame*constraint.v1;
        auto secondAttachment=rigid_structure2->frame*constraint.v2;
        (*vertices)[0].set(firstAttachment(0),firstAttachment(1),firstAttachment(2));
        (*vertices)[1].set(secondAttachment(0),secondAttachment(1),secondAttachment(2));
        lineGeometry->setVertexArray(vertices);
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RELATIVE_POSITION_CONSTRAINT)
GENERIC_TYPE_DEFINITION(RELATIVE_POSITION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(RELATIVE_POSITION_CONSTRAINT,void)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(simulation.data.Find("RIGID_STRUCTURE_DATA"));
    auto relative_position_constraint=simulation.force.template Find_Or_Create<RELATIVE_POSITION_CONSTRAINT<TV>>();
    Json::Value constraints=node["constraints"];
    for(Json::ValueIterator it=constraints.begin();it!=constraints.end();it++){
        typename RELATIVE_POSITION_CONSTRAINT<TV>::CONSTRAINT constraint;
        constraint.s1=rigid_data->Structure_Index((*it)["structure1"].asString());
        Parse_Vector((*it)["offset1"],constraint.v1);
        constraint.s2=rigid_data->Structure_Index((*it)["structure2"].asString());
        Parse_Vector((*it)["offset2"],constraint.v2);
        constraint.target_distance=(*it)["distance"].asDouble();
        relative_position_constraint->constraints.push_back(constraint);
    }
    relative_position_constraint->stored_forces.resize(relative_position_constraint->constraints.size());
    relative_position_constraint->stored_forces.setZero();
    return 0;
}
