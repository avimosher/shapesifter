#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/EQUATION.h>
#include <Equation/TRUST_REGION.h>
#include <Evolution/EVOLUTION.h>
#include <Parsing/PARSE_SCENE.h>
#include <Utilities/LOG.h>
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
    simulation->write=false;
    if(!PARSE_SCENE<TV>::Parse_Scene(std::cin,*simulation)){return 1;}

    auto solver=simulation->evolution.template Find<TRUST_REGION<TV>>();
    auto equation=solver->equation;
    auto& data=simulation->data;
    auto& force=simulation->force;
    T dt=0.1;
    T time=0;
    T epsilon=1;

    //***********************
    // Set up initial state
    //***********************
    Matrix<T,Dynamic,1> positions;
    data.Pack_Positions(positions);
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
        T predicted_delta_f=gradient0.dot(eps*unknowns);//+(T).5*eps*eps*unknowns.transpose()*hessian*unknowns;
        T error=f1-f0-predicted_delta_f;
        return error;
    };

    T last;
    for(int i=0;i<10;i++,epsilon/=2){
        T next=Evaluate_Step_Error(epsilon);
        if(i>1){std::cout<<"Ratio: "<<last/next<<std::endl;}
        last=next;
    }
    return 0;
}
