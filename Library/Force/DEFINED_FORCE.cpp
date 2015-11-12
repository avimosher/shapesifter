#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/DEFINED_FORCE.h>
#include <Force/FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <iostream>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void DEFINED_FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    for(const auto& force : forces){
        right_hand_side.template block<d,1>((t+d)*force.s,0)+=force.f;
        right_hand_side.template block<t,1>((t+d)*force.s+d,0)+=(rigid_data->structures[force.s]->frame.orientation*force.v).cross(force.f);}
    constraint_terms.resize(0,rigid_data->Velocity_DOF());
    constraint_forces.resize(rigid_data->Velocity_DOF(),0);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(DEFINED_FORCE)
DEFINE_AND_REGISTER_PARSER(DEFINED_FORCE,void)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto defined_force=std::make_shared<DEFINED_FORCE<TV>>();
    Json::Value forces=node["forces"];
    for(Json::ValueIterator it=forces.begin();it!=forces.end();it++){
        typename DEFINED_FORCE<TV>::FORCE_POKE force;
        force.s=rigid_data->Structure_Index((*it)["structure"].asString());
        Parse_Vector((*it)["offset"],force.v,TV());
        Parse_Vector((*it)["force"],force.f);
        defined_force->forces.push_back(force);}
    simulation.force.push_back(defined_force);
    return 0;
}
