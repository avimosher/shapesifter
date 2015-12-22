#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/MATRIX_BUNDLE.h>
#include <Force/BROWNIAN_FORCE.h>
#include <Force/FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/LOG.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void BROWNIAN_FORCE<TV>::
Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    SparseMatrix<T>& constraint_forces=system.Matrix_Block(data,force,*rigid_data,*this);
    SparseMatrix<T>& constraint_terms=system.Matrix_Block(data,force,*this,*rigid_data);
    Matrix<T,Dynamic,1>& right_hand_side=system.RHS(data,force,*rigid_data);

    T k_B=1.38e-2; // pN nm/K
    T kT=k_B*temperature; // pN nm
    T eta=data.globals["eta"];
    if(stochastic){
        T one_over_dt=1/dt;
        stored_right_hand_side.resize(right_hand_side.rows(),1);
        stored_right_hand_side.setZero();
        for(int i=0;i<(*rigid_data).structures.size();i++){
            T radius=rigid_data->structures[i]->radius;
            T linear_drag=6*M_PI*eta*radius;
            T translational_diffusion_coefficient=kT/linear_drag;
            T translational_variance=sqrt(2*translational_diffusion_coefficient*dt);
            T random_displacement=data.random.Gaussian(T(),translational_variance);
            TV random_orientation=data.random.template Direction<TV>();
            //stored_right_hand_side.template block<d,1>((t+d)*i,0)=linear_drag*random_displacement*random_orientation*one_over_dt;
            TV linear_force=linear_drag*random_displacement*random_orientation*one_over_dt;
            linear_force[2]=0; // eliminate Z component
            stored_right_hand_side.template block<d,1>((t+d)*i,0)=linear_force;

            LOG::cout<<"Force on "<<rigid_data->structures[i]->name<<": "<<(linear_drag*random_displacement*random_orientation*one_over_dt).transpose()<<" position: "<<rigid_data->structures[i]->frame.position.transpose()<<std::endl;

            /*T rotational_resistance=8*M_PI*eta*std::pow(radius,3);
            T rotational_diffusion_coefficient=kT/rotational_resistance;
            T_SPIN random_spin_orientation=data.random.template Direction<T_SPIN>();
            T rotational_variance=sqrt(2*rotational_diffusion_coefficient*dt);
            T random_angle=data.random.Gaussian(T(),rotational_variance);
            stored_right_hand_side.template block<T_SPIN::SizeAtCompileTime,1>(TWIST<TV>::STATIC_SIZE*i+TV::SizeAtCompileTime,0)=random_spin_orientation*rotational_resistance*random_angle*one_over_dt;*/}}
    right_hand_side+=stored_right_hand_side;
    constraint_terms.resize(0,rigid_data->Velocity_DOF());
    constraint_forces.resize(rigid_data->Velocity_DOF(),0);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(BROWNIAN_FORCE)
DEFINE_AND_REGISTER_PARSER(BROWNIAN_FORCE,void)
{
    auto brownian_force=std::make_shared<BROWNIAN_FORCE<TV>>();
    Parse_Scalar(node["temperature"],brownian_force->temperature);
    simulation.force.push_back(brownian_force);
    return 0;
}
