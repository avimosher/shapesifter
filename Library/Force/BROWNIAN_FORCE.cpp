//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class BROWNIAN_FORCE
///////////////////////////////////////////////////////////////////////
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/BROWNIAN_FORCE.h>
#include <Force/FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <iostream>
#include <math.h>
#include <random>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> BROWNIAN_FORCE<TV>::
BROWNIAN_FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> BROWNIAN_FORCE<TV>::
~BROWNIAN_FORCE()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> void BROWNIAN_FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,std::vector<Triplet<T>>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side)
{
    std::random_device rd; //TODO: separate random class
    std::mt19937 generator(rd());
    T temperature=300; // K
    T one_over_dt=1/dt;
    T k_B=1.38e-2; // pN nm/K
    T kT=k_B*temperature; // pN nm
    T eta=3.5;
    T radius=1;
    for(auto iterator=data.find("RIGID_STRUCTURE_DATA");iterator!=data.end();++iterator){
        auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(iterator->second);
        for(int i=0;i<(*rigid_data).structures.size();i++){
            T translational_diffusion_coefficient=kT/(6*M_PI*eta*radius);
            T translational_variance=sqrt(2*translational_diffusion_coefficient*dt);
            std::normal_distribution<> distribution(0,translational_variance);
            T random_displacement=distribution(generator);
            TV test;test.fill(random_displacement);
std::cout<<test<<std::endl;
            right_hand_side.template block<TV::SizeAtCompileTime,1>(TV::SizeAtCompileTime*i,0)=test;
        }
    }
    //force_terms.push_back(Triplet<T>(0,0,k*dt));
    //right_hand_side(0)=-
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(BROWNIAN_FORCE)
DEFINE_AND_REGISTER_PARSER(BROWNIAN_FORCE)
{
simulation.force.push_back(std::make_shared<BROWNIAN_FORCE<TV>>());
std::cout<<"How unexpected"<<std::endl;
}
