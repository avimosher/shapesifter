#if 0
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/ELECTROSTATIC_FORCE.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void ELECTROSTATIC_FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    std::vector<AlignedBox<T,d>> bounding_list(charges.size());
    for(int s=0;s<charges.size();s++){
        TV charge_location=rigid_data->structures[std::get<0>(charges[s])]->frame*std::get<2>(charges[s]);
        bounding_list[s]=AlignedBox<T,d>(charge_location,charge_location);
        index_list[s]=s;
    }
    KdBVH<T,3,int> hierarchy(index_list.begin(),index_list.end(),bounding_list.begin(),bounding_list.end());
    for(int c1=0;c1<charges.size();c1++){
        int s1=std::get<0>(charges[c1]);
        auto structure1=rigid_data->structures[s1];
        auto charge_location=structure1->frame*std::get<2>(charges[c1]);
        PROXIMITY_SEARCH<TV> proximity_search(data,charge_location,distance_threshold);
        BVIntersect(hierarchy,proximity_search);
        for(auto c2 : proximity_search.candidates){
            int s2=std::get<0>(charges[c2]);
            if(s1==s2){continue;} // don't look at intra-body charge interactions
            // TODO: compute the RHS force and derivative contributions.  Pretty compact.
        }
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ELECTROSTATIC_FORCE<TV>::
Parse_Terms(Json::Value& charges_node,const int body)
{
    for(Json::ValueIterator it=charges_node.begin();it!=charges_node.end();it++){
        CHARGE charge;
        std::get<0>(charge)=body;
        Parse_Scalar((*it)["charge"],std::get<1>(charge));
        Parse_Scalar((*it)["offset"],std::get<2>(charge));
        charges.push_back(charge);
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(ELECTROSTATIC_FORCE)
DEFINE_AND_REGISTER_PARSER(ELECTROSTATIC_FORCE,void)
{
    auto electrostatic_force=simulation.force.template Find_Or_Create<ELECTROSTATIC_FORCE<TV>>();
    // TODO: register parse terms for "charges" with the parser for RIGID_STRUCTURE
}
#endif
