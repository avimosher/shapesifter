#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Equation/MATRIX_BUNDLE.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/RELATIVE_POSITION_CONSTRAINT.h>
#include <Math/Relative_Position_Constraint_Term.h>
#include <Math/Relative_Position_Force.h>
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
Identify_Interactions_And_Compute_Errors(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    VECTOR& right_hand_side=system.RHS(data,force,*rigid_data);
    VECTOR& constraint_right_hand_side=system.RHS(data,force,*this);

    constraint_right_hand_side.resize(constraints.size(),1);
    for(int i=0;i<constraints.size();i++){
        const CONSTRAINT& constraint=constraints[i];
        std::array<int,2> indices={constraint.s1,constraint.s2};
        auto structure1=rigid_data->structures[indices[0]];
        auto structure2=rigid_data->structures[indices[1]];
        std::array<FRAME<TV>,2> frames={structure1->frame,structure2->frame};
        std::array<TV,2> offsets={constraint.v1,constraint.v2};
        std::array<T_SPIN,2> spins={structure1->twist.angular,structure2->twist.angular};
        std::array<TV,2> spun_offsets,base_offsets;
        for(int j=0;j<2;j++){
            spun_offsets[j]=frames[j].orientation*offsets[j];
            base_offsets[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]).inverse()*spun_offsets[j];}

        TV relative_position=data.Minimum_Offset(frames[0]*constraint.v1,frames[1]*constraint.v2);
        Relative_Position_Force<TV>::Right_Hand_Sides(indices,relative_position,stored_forces[i],spins,base_offsets,right_hand_side);
        constraint_right_hand_side[i]=relative_position.norm()-constraint.target_distance;
    }
    // This has to go before Flatten calls - it determines the DOF for this force
    stored_forces.resize(constraints.size());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Compute_Derivatives(DATA<TV>& data,FORCE<TV>& force,MATRIX_BUNDLE<TV>& system)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    std::vector<Triplet<T>>& force_terms=system.Matrix_Block_Terms(data,force,*rigid_data);
    std::vector<Triplet<T>>& hessian_terms=system.Hessian_Block_Terms(data,force,*rigid_data,*rigid_data);
    std::vector<Triplet<T>>& constraint_force_terms=system.Hessian_Block_Terms(data,force,*this,*rigid_data);
    std::vector<Triplet<T>>& force_constraint_terms=system.Hessian_Block_Terms(data,force,*rigid_data,*this);
    const SparseMatrix<T>& f_scaling=system.inverse_inertia_matrices[0];
    VECTOR& constraint_right_hand_side=system.RHS(data,force,*this);
    const Matrix<T,Dynamic,1>& force_balance_error=system.RHS(data,force,*rigid_data);

    std::vector<Triplet<T>> terms;
    std::vector<Triplet<T>> forces;
    for(int i=0;i<constraints.size();i++){
        const CONSTRAINT& constraint=constraints[i];
        std::array<int,2> indices={constraint.s1,constraint.s2};
        auto structure1=rigid_data->structures[indices[0]];
        auto structure2=rigid_data->structures[indices[1]];
        std::array<FRAME<TV>,2> frames={structure1->frame,structure2->frame};
        std::array<T_SPIN,2> spins={structure1->twist.angular,structure2->twist.angular};
        std::array<TV,2> offsets={constraint.v1,constraint.v2};
        std::array<TV,2> spun_offsets,base_offsets;
        for(int j=0;j<2;j++){
            spun_offsets[j]=frames[j].orientation*offsets[j];
            base_offsets[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]).inverse()*spun_offsets[j];}
        TV relative_position=data.Minimum_Offset(frames[0]*constraint.v1,frames[1]*constraint.v2);

        Relative_Position_Constraint_Term<TV>::Constraint_Velocity_Derivatives(indices,i,relative_position,spins,base_offsets,terms);
        Relative_Position_Force<TV>::Balance_Force_Derivatives(indices,i,relative_position,spins,base_offsets,forces);
        Relative_Position_Force<TV>::Balance_Velocity_Derivatives(indices,stored_forces[i],relative_position,spins,base_offsets,force_terms);
        //RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Constraint_Second_Derivatives(force_balance_error,indices,i,constraint_right_hand_side[i],stored_forces[i],relative_position,spins,base_offsets,f_scaling,hessian_terms,force_constraint_terms,constraint_force_terms);
    }
    system.Build_Jacobian_Block(data,force,*this,*rigid_data,terms);
    system.Build_Jacobian_Block(data,force,*rigid_data,*this,forces);
    system.Flatten_Hessian_Block(data,force,*this,*rigid_data,constraint_force_terms);
    system.Flatten_Hessian_Block(data,force,*rigid_data,*this,force_constraint_terms);
}
///////////////////////////////////////////////////////////////////////
#ifdef VIEWER
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
    osg::Group* group=node->asGroup();
    osg::Group* relative_position_group=(osg::Group*)getNamedChild(group,Static_Name());
    if(!relative_position_group){
        relative_position_group=new osg::Group();
        relative_position_group->setName(Static_Name());
        for(int i=0;i<constraints.size();i++){
//            relative_position_group->addChild(createLine(osg::Vec4(1.0f,1.0f,0.0f,1.0f)));}
            relative_position_group->addChild(createLine(osg::Vec4(0.3f,0.3f,0.3f,1.0f)));}
        group->addChild(relative_position_group);
    }
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    for(int i=0;i<constraints.size();i++){
        const CONSTRAINT& constraint=constraints[i];
        auto rigid_structure1=rigid_data->structures[constraint.s1];
        auto rigid_structure2=rigid_data->structures[constraint.s2];
        auto firstAttachment=rigid_structure1->frame*constraint.v1;
        auto secondAttachment=rigid_structure2->frame*constraint.v2;
        updateLine<TV>((osg::Geode*)relative_position_group->getChild(i),{firstAttachment,secondAttachment});
    }
}
#endif
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
