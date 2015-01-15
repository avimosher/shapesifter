#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/VOLUME_EXCLUSION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    T one_over_dt=1/dt;
    RANDOM<TV> random;
    
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.find("RIGID_STRUCTURE_DATA")->second);
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<T> rhs;
    RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
    int constraint_count=0;
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
            //if(constraint_violation<distance_condition || (call_count==remembered_force.x && remembered_force.y<0)){
            std::cout<<"Violation: "<<constraint_violation<<" displacement: "<<distance<<" orientation: "<<direction<<std::endl;
            if(constraint_violation<distance_condition){
                std::cout<<"Try constraint"<<std::endl;
                terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraint_count,s2,direction.transpose()*index_map.Velocity_Map(offset2)));
                terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraint_count,s1,-direction.transpose()*index_map.Velocity_Map(offset1)));
                rhs.push_back(-constraint_violation);
                //constraint_rhs(constraint_count,0)=-constraint_violation;
                constraint_count++;
            }
        }
    }
    constraint_terms.resize(constraint_count,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE*rigid_data->structures.size());
    constraint_rhs.resize(rhs.size(),1);
    for(int i=0;i<rhs.size();i++){constraint_rhs(i,0)=rhs[i];}
    Flatten_Matrix(terms,constraint_terms);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(VOLUME_EXCLUSION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(VOLUME_EXCLUSION_CONSTRAINT,void)
{
    auto volume_exclusion_constraint=std::make_shared<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
    simulation.force.push_back(volume_exclusion_constraint);
    return 0;
}
