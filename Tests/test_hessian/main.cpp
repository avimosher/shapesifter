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
    structure3->frame.position=3*TV::UnitX();
    structure3->name="third";
    structure3->radius=1;
    structure3->Initialize_Inertia(3.5);
    rigid_data->structures.push_back(structure3);

    Matrix<T,Dynamic,1> positions;
    data.Pack_Positions(positions);

    std::cout<<std::setprecision(9);

    auto relative_position_constraint=simulation->force.template Find_Or_Create<RELATIVE_POSITION_CONSTRAINT<TV>>();
    typename RELATIVE_POSITION_CONSTRAINT<TV>::CONSTRAINT constraint;
    constraint.s1=rigid_data->Structure_Index("first");
    constraint.v1.setZero();//TV::UnitY();
    constraint.s2=rigid_data->Structure_Index("second");
    constraint.v2.setZero();//TV::UnitY();
    constraint.target_distance=4;
    relative_position_constraint->constraints.push_back(constraint);

    typename RELATIVE_POSITION_CONSTRAINT<TV>::CONSTRAINT constraint2;
    constraint2.s1=rigid_data->Structure_Index("second");
    constraint2.v1.setZero();//TV::UnitY();
    constraint2.s2=rigid_data->Structure_Index("third");
    constraint2.v2.setZero();//TV::UnitY();
    constraint2.target_distance=4;
    relative_position_constraint->constraints.push_back(constraint2);


    relative_position_constraint->stored_forces.resize(2);
    relative_position_constraint->stored_forces.setZero();


    T dt=0.1;
    T time=0;
    auto equation=new NONLINEAR_EQUATION<TV>();
    equation->Initialize(data,force);
    equation->Linearize(data,force,dt,time,false);
    Matrix<T,Dynamic,1> unknowns=equation->Get_Unknowns(data,force);
    data.random.Direction(unknowns);
    //unknowns.setZero();
    //unknowns(0)=1;
    // block out the rotation terms
    for(int i=0;i<3;i++){
        unknowns[3+i]=0;
        unknowns[9+i]=0;
        unknowns[15+i]=0;
    }
    for(int i=1;i<unknowns.size();i++){unknowns(i)=0;}
    equation->Increment_Unknowns(unknowns,data,force);
    equation->Linearize(data,force,dt,time,false);

    Matrix<T,Dynamic,1> gradient0;equation->Gradient(gradient0);
    T f0=equation->Evaluate();
    SparseMatrix<T> jacobian;
    equation->Jacobian(jacobian);
    T term0=jacobian.coeffRef(0,0)*jacobian.coeffRef(0,0)*unknowns(0);
    std::cout<<"Term 0: "<<term0<<std::endl;

    T epsilon=1e-5;
    TV randdir1,randdir2;
    data.random.Direction(randdir1);
    data.random.Direction(randdir2);
    randdir1*=epsilon;
    randdir2*=epsilon;
    unknowns.setZero();
    for(int i=0;i<3;i++){
        unknowns(i)=randdir1(i);
        unknowns(i+6)=randdir2(i);
    }
    //std::cout<<"Unknowns: "<<std::endl;
    //std::cout<<unknowns.transpose()<<std::endl;
    //unknowns.setZero();unknowns(0)=epsilon;
    //std::cout<<"Term0*delta: "<<term0*epsilon<<std::endl;
    //unknowns(2)=epsilon;
    //unknowns(6)=epsilon;
    //data.random.Direction(unknowns);
    //unknowns*=epsilon;
    // block out the rotation terms
    for(int i=0;i<3;i++){
        unknowns[3+i]=0;
        unknowns[9+i]=0;
        unknowns[15+i]=0;
    }
    std::cout<<unknowns.transpose()<<std::endl;

    
    SparseMatrix<T> hessian;
    equation->Hessian(hessian);
    Matrix<T,Dynamic,1> predicted_delta=hessian*unknowns;
    //std::cout<<hessian<<std::endl;
    //std::cout<<"dx: "<<unknowns.transpose()<<std::endl;


    data.Unpack_Positions(positions);
    equation->Increment_Unknowns(unknowns,data,force);
    equation->Linearize(data,force,dt,time,false);

    Matrix<T,Dynamic,1> gradient1;equation->Gradient(gradient1);
    T f1=equation->Evaluate();
    T predicted_delta_f=gradient0.dot(unknowns);

    std::cout<<"g0: "<<gradient0.transpose()<<std::endl;
    std::cout<<"g1: "<<gradient1.transpose()<<std::endl;
    std::cout<<"f0: "<<f0<<" f1: "<<f1<<" Delta f: "<<f1-f0<<" predicted: "<<predicted_delta_f<<" quality: "<<(f1-f0-predicted_delta_f)/epsilon<<std::endl;

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
