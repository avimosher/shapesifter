#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/RELATIVE_POSITION_CONSTRAINT.h>
#include <Utilities/RANDOM.h>
#include <algorithm>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <string>

#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/NONLINEAR_EQUATION.h>
using namespace Mechanics;

int main(int argc,char **argv)
{
    typedef double T;
    typedef Matrix<T,3,1> TV;
    auto simulation=std::make_shared<SIMULATION<TV>>();
    DATA<TV>& data=simulation->data;
    FORCE<TV>& force=simulation->force;

    auto rigid_data=data.template Find_Or_Create<RIGID_STRUCTURE_DATA<TV>>();
    auto structure1=std::make_shared<RIGID_STRUCTURE<TV>>();
    structure1->frame.position=TV();
    structure1->name="first";
    structure1->radius=1;
    structure1->Initialize_Inertia(3.5);
    rigid_data->structures.push_back(structure1);
    auto structure2=std::make_shared<RIGID_STRUCTURE<TV>>();
    structure2->frame.position=2*TV::UnitX();
    structure2->name="second";
    structure2->radius=1;
    structure2->Initialize_Inertia(3.5);
    structure2->twist.linear[0]=0;
    rigid_data->structures.push_back(structure2);


    auto relative_position_constraint=simulation->force.template Find_Or_Create<RELATIVE_POSITION_CONSTRAINT<TV>>();
    typename RELATIVE_POSITION_CONSTRAINT<TV>::CONSTRAINT constraint;
    constraint.s1=rigid_data->Structure_Index("first");
    constraint.v1=TV::UnitY();
    constraint.s2=rigid_data->Structure_Index("second");
    constraint.v2=TV::UnitY();
    constraint.target_distance=2;
    relative_position_constraint->constraints.push_back(constraint);
    relative_position_constraint->stored_forces.resize(1);
    relative_position_constraint->stored_forces.setZero();


    T dt=0.1;
    T time=0;
    auto equation=new NONLINEAR_EQUATION<TV>();
    equation->Initialize(data,force);
    equation->Linearize(data,force,dt,time,false);
    Matrix<T,Dynamic,1> unknowns=equation->Get_Unknowns(data,force);
    data.random.Direction(unknowns);
    equation->Increment_Unknowns(unknowns,data,force);
    equation->Linearize(data,force,dt,time,false);

    Matrix<T,Dynamic,1> gradient0;equation->Gradient(gradient0);

    T epsilon=1e-6;
    data.random.Direction(unknowns);
    unknowns*=epsilon;

    
    SparseMatrix<T> hessian;
    equation->Hessian(hessian);
    Matrix<T,Dynamic,1> predicted_delta=hessian*unknowns;


    equation->Increment_Unknowns(unknowns,data,force);
    equation->Linearize(data,force,dt,time,false);

    Matrix<T,Dynamic,1> gradient1;equation->Gradient(gradient1);

    std::cout<<(gradient1-gradient0).transpose()<<std::endl;
    std::cout<<predicted_delta.transpose()<<std::endl;
    std::cout<<"Quality: "<<(gradient1-gradient0-predicted_delta).norm()/epsilon<<std::endl;

    /*
    Matrix<T,Dynamic,1> positions;
    simulation->data.Pack_Positions(positions);

    // compare to actual results
    Matrix<T,Dynamic,1> initial_velocities;initial_velocities.setZero();
    structure2->twist.angular[2]=.1;
    simulation->data.Pack_Velocities(initial_velocities);
    simulation->data.Step();

    SparseMatrix<T> gradient;
    relative_position_constraint->Special(simulation->data,1,1,gradient);


    // get actuall new constraint value
    Matrix<T,Dynamic,1> satisfaction;
    relative_position_constraint->Constraint_Satisfaction(simulation->data,1,1,satisfaction);
    std::cout<<"Satisfaction: "<<satisfaction<<std::endl;

    simulation->data.Unpack_Positions(positions);
    RANDOM<T> random;
    structure2->twist.angular+=random.template Direction<TV>()/10;
    Matrix<T,Dynamic,1> final_velocities;final_velocities.setZero();
    simulation->data.Pack_Velocities(final_velocities);


    std::cout<<"Gradient: "<<gradient<<std::endl;
    //std::cout<<"Velocities: "<<velocities.transpose()<<std::endl;
    std::cout<<"Predicted: "<<satisfaction+gradient*(final_velocities-initial_velocities)<<std::endl;
    simulation->data.Step();
    relative_position_constraint->Constraint_Satisfaction(simulation->data,1,1,satisfaction);
    std::cout<<"Satisfaction: "<<satisfaction<<std::endl;
    */
    return 0;
}
