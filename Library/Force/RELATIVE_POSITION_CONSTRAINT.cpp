#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/RELATIVE_POSITION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/LOG.h>
#include <Utilities/MATH.h>
#include <Utilities/OSG_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
#include <osg/Geometry>
#include <osg/Geode>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
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
        auto structure1=rigid_data->structures[body_index1];
        auto structure2=rigid_data->structures[body_index2];
        FRAME<TV> frame1=structure1->frame;
        FRAME<TV> frame2=structure2->frame;
        TV direction=data.Minimum_Offset(frame1*constraint.v1,frame2*constraint.v2);
        T distance=direction.norm();
        /*LOG::cout<<"S1 position: "<<frame1.position.transpose()<<" constraint point: "<<(frame1*constraint.v1).transpose()<<std::endl;
          LOG::cout<<"S2 position: "<<frame2.position.transpose()<<" constraint point: "<<(frame2*constraint.v2).transpose()<<std::endl;*/
        //LOG::cout<<"Constraint distance: "<<distance<<" between "<<structure1->name<<" and "<<structure2->name<<" direction "<<direction.transpose()<<" force "<<stored_forces(i)<<std::endl;
        LOG::cout<<"Constraint distance: "<<distance<<" between "<<body_index1<<" and "<<body_index2<<" direction "<<direction.transpose()<<" force "<<stored_forces(i)<<" v1: "<<constraint.v1.transpose()<<" v2: "<<constraint.v2.transpose()<<std::endl;
        TV x1=frame1*constraint.v1;
        TV x2=frame2*constraint.v2;
        CONSTRAINT_VECTOR DC_DA2=RIGID_STRUCTURE_INDEX_MAP<TV>::DC_DA(*structure2,constraint.v2,x1,x2,direction);
        CONSTRAINT_VECTOR DC_DA1=RIGID_STRUCTURE_INDEX_MAP<TV>::DC_DA(*structure1,constraint.v1,x1,x2,direction);
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index2,DC_DA2));
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index1,-DC_DA1));
        Matrix<T,t+d,t+d> force_balance_contribution2=stored_forces(i)*RIGID_STRUCTURE_INDEX_MAP<TV>::DF_DA(*structure2,constraint.v2,x1,x2,direction);
        Matrix<T,t+d,t+d> force_balance_contribution1=stored_forces(i)*RIGID_STRUCTURE_INDEX_MAP<TV>::DF_DA(*structure1,constraint.v1,x1,x2,direction); // the two negatives
                                                                                                                                                               // actually cancel out.  No
                                                                                                                                                               // net negative sign
        TV offset1=frame1.orientation*constraint.v1;
        TV offset2=frame2.orientation*constraint.v2;
        constraint_rhs[i]=(constraint.target_distance-distance);
        //LOG::cout<<"F_i: "<<constraint_rhs[i]<<std::endl;
#if 1
        for(int j=0;j<t+d;j++){
            for(int k=0;k<t+d;k++){
                if(fabs(force_balance_contribution1(j,k))>1e-6){
                    force_terms.push_back(Triplet<T>(body_index1*(t+d)+j,body_index1*(t+d)+k,force_balance_contribution1(j,k)));
                }
                if(fabs(force_balance_contribution2(j,k))>1e-6){
                    force_terms.push_back(Triplet<T>(body_index2*(t+d)+j,body_index2*(t+d)+k,force_balance_contribution2(j,k)));
                }
            }
        }
#endif
        // contribution to force-balance RHS
        //LOG::cout<<"REL RHS contribution: "<<(DC_DA1.transpose()*stored_forces[i]).transpose()<<std::endl;
        right_hand_side.template block<t+d,1>(body_index1*(t+d),0)+=DC_DA1.transpose()*stored_forces[i];
        right_hand_side.template block<t+d,1>(body_index2*(t+d),0)-=DC_DA2.transpose()*stored_forces[i];
    }
    constraint_terms.resize(constraints.size(),RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE*rigid_data->structures.size());
    Flatten_Matrix(terms,constraint_terms);
    stored_forces.resize(constraints.size());
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
