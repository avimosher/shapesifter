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
    T rotational_drag=8*M_PI*eta*std::pow(radius,3);
    T linear_drag=6*M_PI*eta*radius;
    if(!collision_extent){
        moi.translation=TV::Constant(linear_drag);
        moi.rotation=T_SPIN::Constant(rotational_drag);
    }
    else{
        T p=collision_extent/collision_radius+1;
        T E=sqrt(abs(p*p-1))/p;
        T S=2*atanh(E)/E;
        
        // Perrin friction factors
        T half_length=collision_extent+collision_radius;
        T equivalent_linear_drag=6*M_PI*eta*std::pow(half_length*collision_radius*collision_radius,1.0f/3.0f);
        moi.translation=TV::Constant(2*std::pow(p,2.0f/3.0f)/S*equivalent_linear_drag);
        moi.rotation(2)=4.0f/3.0f*E*E/(2-S/(p*p))*rotational_drag;
        for(int spin=0;spin<2;spin++){
            moi.rotation(spin)=4.0f/3.0f*(1/(p*p)-p*p)/(2-S*(2-1/(p*p)))*rotational_drag;
        }
    }
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RIGID_STRUCTURE)
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE,void)
{
    std::shared_ptr<RIGID_STRUCTURE<TV>> structure=std::make_shared<RIGID_STRUCTURE<TV>>();
    Parse_Vector(node["position"],structure->frame.position,TV());
    Parse_Scalar(node["collision_extent"],structure->collision_extent);
    Parse_Scalar(node["radius"],structure->radius,(T)0);
    Parse_Scalar(node["collision_radius"],structure->collision_radius,structure->radius);
    structure->name=node["name"].asString();
    auto data_element=simulation.data.Find(RIGID_STRUCTURE_DATA<TV>::Static_Name());
    auto rigid_structure_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data_element);
    rigid_structure_data->structures.push_back(structure);
    return 0;
}
