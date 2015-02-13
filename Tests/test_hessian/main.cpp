#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/RELATIVE_POSITION_CONSTRAINT.h>
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
    structure2->frame.position=TV::UnitX();
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

    relative_position_constraint->Special(simulation->data,1,1);
    return 0;
}
