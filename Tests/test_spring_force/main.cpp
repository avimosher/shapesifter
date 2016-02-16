#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/SPRING_FORCE.h>
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

    Matrix<T,Dynamic,1> positions;
    data.Pack_Positions(positions);

    auto spring_force=simulation->force.template Find_Or_Create<SPRING_FORCE<TV>>();
    typename SPRING_FORCE<TV>::SPRING spring;
    spring.index={0,1};
    spring.offset={TV::UnitX(),TV::UnitY()};
    //spring.offset={TV(),TV()};
    spring.target_distance=3;
    spring.stiffness=10;
    spring_force->springs.push_back(spring);

    T dt=0.1;
    T time=0;
    T epsilon=1e-1;
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
    /*for(int i=0;i<3;i++){
        unknowns[3+i]=0;
        unknowns[9+i]=0;
    }*/

    auto Evaluate_Step_Error = [&](T eps){

        data.Unpack_Positions(positions);
        equation->Increment_Unknowns(eps*unknowns,data,force);
        equation->Linearize(data,force,dt,time,false);
        equation->Increment_Unknowns(-eps*unknowns,data,force);

        T f1=equation->Evaluate();
        T predicted_delta_f=gradient0.dot(eps*unknowns);//+(T).5*eps*eps*unknowns.transpose()*hessian*unknowns;
        T error=f1-f0-predicted_delta_f;
        return error;
    };

    T last;
    for(int i=0;i<6;i++,epsilon/=2){
        T next=Evaluate_Step_Error(epsilon);
        if(i>1){
            std::cout<<"Ratio: "<<last/next<<std::endl;
            //std::cout<<"Ratio: "<<Evaluate_Step_Error(epsilon)/Evaluate_Step_Error(epsilon/2)<<std::endl;
        }
        last=next;
    }
    return 0;
}
