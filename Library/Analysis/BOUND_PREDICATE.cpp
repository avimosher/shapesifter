#include <Analysis/BOUND_PREDICATE.h>
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
#include <Force/FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar BOUND_PREDICATE<TV>::
Scalar(const SIMULATION<TV>& simulation)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto association_dissociation_constraint=simulation.force.template Find_Or_Create<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>();
    int first_binder_index=rigid_data->Structure_Index(first_binder);
    int second_binder_index=rigid_data->Structure_Index(second_binder);
    for(auto& constraint : association_dissociation_constraint->constraints){
        auto interaction_type=association_dissociation_constraint->interaction_types[std::get<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::CONSTRAINT_INTERACTION>(constraint)];
        auto first_site=interaction_type.sites[0][std::get<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::CONSTRAINT_BODY1>(constraint)];
        auto second_site=interaction_type.sites.back()[std::get<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::CONSTRAINT_BODY2>(constraint)];
        int body1=std::get<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::SITE_INDEX>(first_site);
        int body2=std::get<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::SITE_INDEX>(second_site);
        //std::cout<<"body1: "<<body1<<" body2: "<<body2<<" i1: "<<first_binder_index<<" i2: "<<second_binder_index<<std::endl;
        if(std::min(body1,body2)==std::min(first_binder_index,second_binder_index) && std::max(body1,body2)==std::max(first_binder_index,second_binder_index)){
            return 1;
        }
    }
    return 0;
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(BOUND_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(BOUND_PREDICATE,PREDICATE)
{
    auto predicate=std::make_shared<BOUND_PREDICATE<TV>>();
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    Parse_String(node["first_binder"],predicate->first_binder);
    Parse_String(node["second_binder"],predicate->second_binder);
    return predicate;
}
