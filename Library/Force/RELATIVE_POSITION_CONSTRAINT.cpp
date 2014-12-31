#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/RELATIVE_POSITION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
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
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side)
{
    T one_over_dt=1/dt;
    RANDOM<TV> random;
    for(auto iterator=data.find("RIGID_STRUCTURE_DATA");iterator!=data.end();++iterator){
        auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(iterator->second);
        typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
        std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
        RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
        typedef int CONSTRAINT;
        std::vector<CONSTRAINT> constraints;
        for(int i=0;i<constraints.size();i++){
            //const CONSTRAINT& constraint=constraints(i);
            FRAME<TV> frame1;//=rigid_data->Updated_Frame(rigid_structure1->frame,rigid_velocity->twist(body_index1));
            FRAME<TV> frame2;//=rigid_data->Updated_Frame(rigid_structure2->frame,rigid_velocity->twist(body_index2));
            int body_index1;
            int body_index2;
            TV offset1;//=frame1.rotation._transformVector(constraint.structure1.y);
            TV offset2;//=frame2.rotation._transformVector(constraint.structure2.y);
            TV direction;//TODO=data.Minimum_Offset(frame1.position+offset1,frame2.position+offset2);
            T distance=direction.norm();
            direction.normalize();
            terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index2,-direction.transpose()*index_map.Velocity_Map(offset2)));
            terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index1,direction.transpose()*index_map.Velocity_Map(offset1)));
            // right_hand_side(i)=F0; // need constraint RHS
        }
        // TODO: define blocks-per-force
        //matrix.setFromTriplets(terms);
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RELATIVE_POSITION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(RELATIVE_POSITION_CONSTRAINT)
{
    auto relative_position_constraint=std::make_shared<RELATIVE_POSITION_CONSTRAINT<TV>>();
    simulation.force.push_back(relative_position_constraint);
}
