#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
#include <Utilities/RANDOM.h>
#include <algorithm>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <iomanip>
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
    rigid_data->structures.push_back(structure2);

    auto structure3=std::make_shared<RIGID_STRUCTURE<TV>>();
    structure3->frame.position=TV()+3*TV::UnitY();
    structure3->name="third";
    structure3->radius=1;
    structure3->Initialize_Inertia(3.5);
    rigid_data->structures.push_back(structure3);

    auto structure4=std::make_shared<RIGID_STRUCTURE<TV>>();
    structure4->frame.position=2*TV::UnitX()+3*TV::UnitY();
    structure4->name="fourth";
    structure4->radius=1;
    structure4->Initialize_Inertia(3.5);
    rigid_data->structures.push_back(structure4);

    Matrix<T,Dynamic,1> positions;
    data.Pack_Positions(positions);

    auto association_dissociation_constraint=simulation->force.template Find_Or_Create<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>();
    typename ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::INTERACTION_TYPE interaction;
    interaction.bond_distance_threshold=1.5;
    interaction.bond_orientation_threshold=1.5;
    interaction.base_dissociation_time=200;
    interaction.base_association_time=0.001;
    interaction.sites.resize(2);
    interaction.site_offsets.resize(2);
    interaction.sites[0].push_back(std::pair<int,bool>(0,false));
    interaction.site_offsets[0]=TV::UnitX();
    interaction.sites[1].push_back(std::pair<int,bool>(1,false));
    interaction.site_offsets[1]=-TV::UnitX();
    interaction.relative_orientation=ROTATION<TV>();
    association_dissociation_constraint->interaction_types.push_back(interaction);

    T dt=0.1;
    T time=0;
    T epsilon=4e-3;
    auto equation=new NONLINEAR_EQUATION<TV>();

    //***********************
    // Set up initial state
    //***********************

    equation->Initialize(data,force);
    equation->Linearize(data,force,dt,time,true);
    Matrix<T,Dynamic,1> unknowns=equation->Get_Unknowns(data,force);
    data.random.Direction(unknowns);
    equation->Increment_Unknowns(unknowns,data,force);
    equation->Linearize(data,force,dt,time,false);

    Matrix<T,Dynamic,1> gradient0;equation->Gradient(gradient0);
    T f0=equation->Evaluate();
    SparseMatrix<T> hessian;equation->Hessian(hessian);

    //***********************
    // Calculate step result
    //***********************
    data.random.Direction(unknowns);

    auto Evaluate_Step_Error = [&](T eps){

        data.Unpack_Positions(positions);
        equation->Increment_Unknowns(eps*unknowns,data,force);
        equation->Linearize(data,force,dt,time,false);
        equation->Increment_Unknowns(-eps*unknowns,data,force);

        T f1=equation->Evaluate();
        T predicted_delta_f=gradient0.dot(eps*unknowns)+(T).5*eps*eps*unknowns.transpose()*hessian*unknowns;
        T error=f1-f0-predicted_delta_f;
        std::cout<<"Error: "<<error<<std::endl;
        return error;
    };

    std::cout<<"Ratio: "<<Evaluate_Step_Error(epsilon)/Evaluate_Step_Error(epsilon/2)<<std::endl;
    return 0;
}
