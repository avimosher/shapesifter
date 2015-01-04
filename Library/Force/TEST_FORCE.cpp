//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class TEST_FORCE
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/TEST_FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <iostream>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> TEST_FORCE<TV>::
TEST_FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> TEST_FORCE<TV>::
~TEST_FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TEST_FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    // let's say the test force is F=-kv, and v is our only system variable (?)
    // so v=v*-k*dt*v
    for(auto& test_data : data) {
        T k=4;
        std::cout<<"blah"<<std::endl;
        force_terms.push_back(Triplet<T>(0,0,k*dt));
        //right_hand_side(0)=-
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(TEST_FORCE)
DEFINE_AND_REGISTER_PARSER(TEST_FORCE)
{
    simulation.force.push_back(std::make_shared<TEST_FORCE<TV>>());
}
