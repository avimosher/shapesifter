#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE<TV>::
Initialize_Inertia(const T eta)
{
    T radius=substructures[0].radius;
    const TV& capsule_extent=substructures[0].capsule_extent;
    T rotational_drag=8*M_PI*eta*std::pow(radius,3);
    T linear_drag=6*M_PI*eta*radius;
    if(!capsule_extent.norm()){
#if 0
        moi.translation=TV::Constant(1);
        moi.rotation=T_SPIN::Constant(1);
#else
        moi.translation=TV::Constant(linear_drag);
        moi.rotation=T_SPIN::Constant(rotational_drag);
#endif
    }
    else{
        T p=capsule_extent.norm()/radius+1;
        T E=sqrt(fabs(p*p-1))/p;
        T S=2*atanh(E)/E;
        
        // Perrin friction factors
        T half_length=capsule_extent.norm()+radius;
        T equivalent_linear_drag=6*M_PI*eta*std::pow(half_length*radius*radius,1.0f/3.0f);
        moi.translation=TV::Constant(2*std::pow(p,2.0f/3.0f)/S*equivalent_linear_drag);
        moi.rotation(2)=4.0f/3.0f*E*E/(2-S/(p*p))*rotational_drag;
        for(int spin=0;spin<2;spin++){
            moi.rotation(spin)=4.0f/3.0f*(1/(p*p)-p*p)/(2-S*(2-1/(p*p)))*rotational_drag;}}
    initialized=true;
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RIGID_STRUCTURE)
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE,void)
{
    auto structure=std::make_shared<RIGID_STRUCTURE<TV>>();
    Parse_Vector(node["position"],structure->frame.position,TV());
    if(node["substructures"].isNull()){
        structure->substructures.resize(1);
        structure->Substructure(0).Parse(node);}
    else{
        Json::Value substructures=node["substructures"];
        int count=substructures.size();
        structure->substructures.resize(count);
        for(int i=0;i<count;i++){structure->substructures[i].Parse(substructures[i]);}}
    Parse_Scalar(node["kinematic"],structure->kinematic,false);
    Parse_Rotation(node["orientation"],structure->frame.orientation);
    structure->name=node["name"].asString();
    structure->Initialize_Inertia(simulation.data.globals["eta"]);
    auto rigid_structure_data=simulation.data.template Find_Or_Create<RIGID_STRUCTURE_DATA<TV>>();
    rigid_structure_data->structures.push_back(structure);
    LOG::cout<<"Name: "<<structure->name<<" index: "<<rigid_structure_data->structures.size()-1<<std::endl;
    return 0;
}
