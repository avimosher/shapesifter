#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/BROWNIAN_FORCE.h>
#include <Force/FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void BROWNIAN_FORCE<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    T temperature=300; // K
    T one_over_dt=1/dt;
    T k_B=1.38e-2; // pN nm/K
    T kT=k_B*temperature; // pN nm
    T eta=3.5;
    T radius=1;
    RANDOM<TV> random;
    if(stochastic){
        stored_right_hand_side.resize(right_hand_side.rows(),1);
        stored_right_hand_side.setZero();
        for(auto iterator=data.find("RIGID_STRUCTURE_DATA");iterator!=data.end();++iterator){
            auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(iterator->second);
            for(int i=0;i<(*rigid_data).structures.size();i++){
                T translational_diffusion_coefficient=kT/(6*M_PI*eta*radius);
                T translational_variance=sqrt(2*translational_diffusion_coefficient*dt);
                //std::normal_distribution<> distribution(0,translational_variance);
                //T random_displacement=distribution(generator);
                //TV test;test.fill(random_displacement);
                //std::cout<<test<<std::endl;
                stored_right_hand_side.template block<TV::SizeAtCompileTime,1>(TWIST<TV>::STATIC_SIZE*i,0)=random.Direction();
            }
        }
    }
    right_hand_side+=stored_right_hand_side;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(BROWNIAN_FORCE)
DEFINE_AND_REGISTER_PARSER(BROWNIAN_FORCE,void)
{
    auto brownian_force=std::make_shared<BROWNIAN_FORCE<TV>>();
    brownian_force->temperature=node["temperature"].asDouble();
    simulation.force.push_back(brownian_force);
    return 0;
}
