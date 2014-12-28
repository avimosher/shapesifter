#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/RELATIVE_POSITION_CONSTRAINT.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> RELATIVE_POSITION_CONSTRAINT<TV>::
RELATIVE_POSITION_CONSTRAINT()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> RELATIVE_POSITION_CONSTRAINT<TV>::
~RELATIVE_POSITION_CONSTRAINT()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RELATIVE_POSITION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side)
{
    T one_over_dt=1/dt;
    RANDOM<TV> random;
    for(auto iterator=data.find("RIGID_STRUCTURE_DATA");iterator!=data.end();++iterator){
        auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(iterator->second);
        std::vector<Triplet<TRANSLATION_MATRIX>> translation;
        std::vector<Triplet<ROTATION_MATRIX>> rotation;
        Eigen::DiagonalMatrix<TV_TRANSPOSE,Dynamic> projection;
        for(int i=0;i<constraints.size();i++){
            const CONSTRAINT& constraint=constraints(i);
            FRAME<TV> frame1=data.Updated_Frame(rigid_structure1->frame,rigid_velocity->twist(body_index1));
            FRAME<TV> frame2=data.Updated_Frame(rigid_structure2->frame,rigid_velocity->twist(body_index2));
            TV offset1=frame1.rotation._transformVector(constraint.structure1.y);
            TV offset2=frame2.rotation._transformVector(constraint.structure2.y);
            TV direction=data.Minimum_Offset(frame1.position+offset1,frame2.position+offset2);
            T distance=direction.Normalize();
            index_map.Relative_Velocity_Map(body_index1,offset1,body_index2,offset2,translation,rotation);
            for(int axis=0;axis<TV::RowsAtCompileTime;axis++){projection(i,i)(0,axis)=direction(axis);}
            // right_hand_side(i)=F0; // need constraint RHS
        }
        projection*translation;projection*rotation;
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RELATIVE_POSITION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(RELATIVE_POSITION_CONSTRAINT)
{
    auto relative_position_constraint=std::make_shared<RELATIVE_POSITION_CONSTRAINT<TV>>();
    simulation.force.push_back(relative_position_constraint);
}
