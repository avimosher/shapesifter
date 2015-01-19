#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/VOLUME_EXCLUSION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/HASHING.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Unpack_Forces(const Matrix<T,Dynamic,1>& forces)
{
    for(int i=0;i<constraints.size();i++){
        force_memory[constraints[i]]=std::pair<int,T>(call_count,forces(i,0));
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.find("RIGID_STRUCTURE_DATA")->second);
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<T> rhs;
    RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
    constraints.clear();
    for(int s1=0;s1<rigid_data->structures.size();s1++){
        for(int s2=s1+1;s2<rigid_data->structures.size();s2++){
            auto rigid_structure1=rigid_data->structures[s1];
            auto rigid_structure2=rigid_data->structures[s2];
            TV offset1,offset2;
            TV direction=rigid_structure1->Displacement(data,rigid_structure2,offset1,offset2);
            T distance=direction.norm();
            direction.normalize();
            T constraint_violation=distance-rigid_structure1->collision_radius-rigid_structure2->collision_radius;
            T distance_condition=-.001;
            std::pair<int,T> remembered=force_memory[CONSTRAINT(s1,s2)];
            if(constraint_violation<distance_condition || (remembered.first==call_count && remembered.second<0)){
                terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s2,direction.transpose()*index_map.Velocity_Map(offset2)));
                terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s1,-direction.transpose()*index_map.Velocity_Map(offset1)));
                rhs.push_back(-constraint_violation);
                constraints.push_back(CONSTRAINT(s1,s2));
            }
        }
    }
    stored_forces.resize(constraints.size(),1);
    for(int i=0;i<constraints.size();i++){
        stored_forces(i,0)=force_memory[constraints[i]].second;
    }
    constraint_terms.resize(constraints.size(),RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE*rigid_data->structures.size());
    constraint_rhs.resize(rhs.size(),1);
    for(int i=0;i<rhs.size();i++){constraint_rhs(i,0)=rhs[i];}
    Flatten_Matrix(terms,constraint_terms);
    call_count++;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(VOLUME_EXCLUSION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(VOLUME_EXCLUSION_CONSTRAINT,void)
{
    auto volume_exclusion_constraint=std::make_shared<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
    simulation.force.push_back(volume_exclusion_constraint);
    return 0;
}
