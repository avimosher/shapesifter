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
#include <iostream>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    typedef Matrix<T,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE,1> FORCE_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<Triplet<FORCE_VECTOR>> forces;
    constraint_rhs.resize(constraints.size(),1);
    for(int i=0;i<constraints.size();i++){
        const CONSTRAINT& constraint=constraints[i];
        std::vector<int> indices={constraint.s1,constraint.s2};
        auto structure1=rigid_data->structures[indices[0]];
        auto structure2=rigid_data->structures[indices[1]];
        std::vector<FRAME<TV>> frames={structure1->frame,structure2->frame};
        std::vector<T_SPIN> spins={structure1->twist.angular,structure2->twist.angular};
        std::vector<TV> rotated_offsets={frames[0].orientation*constraint.v1,frames[1].orientation*constraint.v2};
        TV relative_position=data.Minimum_Offset(frames[0]*constraint.v1,frames[1]*constraint.v2);

        for(int s=0,sgn=-1;s<2;s++,sgn+=2){
            terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,indices[s],sgn*RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(spins[s],rotated_offsets[s],relative_position)));
            // contribution to force-balance RHS
            FORCE_VECTOR force_direction=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(rotated_offsets[s]).transpose()*relative_position.normalized(); // TODO: may be problematic for distance=0
            right_hand_side.template block<t+d,1>(indices[s]*(t+d),0)-=sgn*force_direction*stored_forces[i];
            forces.push_back(Triplet<FORCE_VECTOR>(indices[s],i,sgn*force_direction));}

        constraint_rhs[i]=(constraint.target_distance-relative_position.norm());
        RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Constraint_Force_Derivatives(indices,stored_forces[i],relative_position,rotated_offsets,spins,force_terms);}
    constraint_terms.resize(constraints.size(),rigid_data->Velocity_DOF());
    Flatten_Matrix(terms,constraint_terms);
    constraint_forces.resize(rigid_data->Velocity_DOF(),constraints.size());
    Flatten_Matrix(forces,constraint_forces);
    stored_forces.resize(constraints.size());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
    osg::Group* group=node->asGroup();
    osg::Group* relative_position_group=(osg::Group*)getNamedChild(group,Static_Name());
    osg::Group* relative_position_error_group=(osg::Group*)getNamedChild(group,"relative_position_error_group");
    if(!relative_position_group){
        relative_position_group=new osg::Group();
        relative_position_group->setName(Static_Name());
        relative_position_error_group=new osg::Group();
        relative_position_error_group->setName("relative_position_error_group");
        for(int i=0;i<constraints.size();i++){
            relative_position_group->addChild(createLine(osg::Vec4(1.0f,1.0f,0.0f,1.0f)));
            relative_position_error_group->addChild(createLine(osg::Vec4(1.0f,1.0f,1.0f,1.0f)));}
        group->addChild(relative_position_group);
        group->addChild(relative_position_error_group);
    }
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    for(int i=0;i<constraints.size();i++){
        const CONSTRAINT& constraint=constraints[i];
        int body_index1=constraint.s1;
        int body_index2=constraint.s2;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        auto firstAttachment=rigid_structure1->frame*constraint.v1;
        auto secondAttachment=rigid_structure2->frame*constraint.v2;
        std::vector<TV> points={firstAttachment,secondAttachment};
        updateLine((osg::Geode*)relative_position_group->getChild(i),points);

        if(errors.rows()==constraints.size()){
            auto meanAttachment=(firstAttachment+secondAttachment)/2;
            auto offsetAttachment=meanAttachment+TV::Unit(2)*errors(i);
            std::vector<TV> error_points={meanAttachment,offsetAttachment};
            updateLine((osg::Geode*)relative_position_error_group->getChild(i),error_points);}
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RELATIVE_POSITION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(RELATIVE_POSITION_CONSTRAINT,void)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
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
