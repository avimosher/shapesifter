#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/TEST_FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void TEST_FORCE<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    // let's say the test force is F=-kv, and v is our only system variable (?)
    // so v=v*-k*dt*v
    for(auto& test_data : data) {
        T k=4;
        //force_terms.push_back(Triplet<T>(0,0,k*dt));
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(TEST_FORCE)
DEFINE_AND_REGISTER_PARSER(TEST_FORCE,void)
{
    simulation.force.push_back(std::make_shared<TEST_FORCE<TV>>());
    return 0;
}
