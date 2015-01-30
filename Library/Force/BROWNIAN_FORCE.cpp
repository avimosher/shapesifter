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
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    T one_over_dt=1/dt;
    T k_B=1.38e-2; // pN nm/K
    T kT=k_B*temperature; // pN nm
    RANDOM<T> random;
    if(stochastic){
        T one_over_dt=1/dt;
        stored_right_hand_side.resize(right_hand_side.rows(),1);
        stored_right_hand_side.setZero();
        auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.Find("RIGID_STRUCTURE_DATA"));
        for(int i=0;i<(*rigid_data).structures.size();i++){
            T radius=rigid_data->structures[i]->radius;
            T linear_drag=6*M_PI*eta*radius;
            T translational_diffusion_coefficient=kT/linear_drag;
            T translational_variance=sqrt(2*translational_diffusion_coefficient*dt);
            T random_displacement=random.Gaussian(T(),translational_variance);
            TV random_orientation=random.template Direction<TV>();
            stored_right_hand_side.template block<TV::SizeAtCompileTime,1>(TWIST<TV>::STATIC_SIZE*i,0)=linear_drag*random_displacement*random_orientation*one_over_dt;
            
            T rotational_resistance=8*M_PI*eta*std::pow(radius,3);
            T rotational_diffusion_coefficient=kT/rotational_resistance;
            T_SPIN random_spin_orientation=random.template Direction<T_SPIN>();
            T rotational_variance=sqrt(2*rotational_diffusion_coefficient*dt);
            T random_angle=random.Gaussian(T(),rotational_variance);
            stored_right_hand_side.template block<T_SPIN::SizeAtCompileTime,1>(TWIST<TV>::STATIC_SIZE*i+TV::SizeAtCompileTime,0)=random_spin_orientation*rotational_resistance*random_angle*one_over_dt;
        }
        std::cout<<stored_right_hand_side.transpose()<<std::endl;
    }
    right_hand_side+=stored_right_hand_side;
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(BROWNIAN_FORCE)
GENERIC_TYPE_DEFINITION(BROWNIAN_FORCE)
DEFINE_AND_REGISTER_PARSER(BROWNIAN_FORCE,void)
{
    auto brownian_force=std::make_shared<BROWNIAN_FORCE<TV>>();
    Parse_Scalar(node["temperature"],brownian_force->temperature);
    brownian_force->eta=simulation.data.globals["eta"];
    simulation.force.push_back(brownian_force);
    return 0;
}
