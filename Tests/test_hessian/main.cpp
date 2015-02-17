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
using namespace Mechanics;

int main(int argc,char **argv)
{
    typedef double T;
    typedef Matrix<T,3,1> TV;
    auto simulation=std::make_shared<SIMULATION<TV>>();

    auto rigid_data=std::make_shared<RIGID_STRUCTURE_DATA<TV>>(); 
    simulation->data.push_back(rigid_data);
    auto structure1=std::make_shared<RIGID_STRUCTURE<TV>>();
    structure1->frame.position=TV();
    structure1->name="first";
    rigid_data->structures.push_back(structure1);
    auto structure2=std::make_shared<RIGID_STRUCTURE<TV>>();
    structure2->frame.position=2*TV::UnitX();
    structure2->name="second";
    rigid_data->structures.push_back(structure2);

    auto relative_position_constraint=simulation->force.template Find_Or_Create<RELATIVE_POSITION_CONSTRAINT<TV>>();
    typename RELATIVE_POSITION_CONSTRAINT<TV>::CONSTRAINT constraint;
    constraint.s1=rigid_data->Structure_Index("first");
    constraint.v1=TV::UnitY();
    constraint.s2=rigid_data->Structure_Index("second");
    constraint.v2=TV::UnitY();
    constraint.target_distance=2;
    relative_position_constraint->constraints.push_back(constraint);

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

    return 0;
}
