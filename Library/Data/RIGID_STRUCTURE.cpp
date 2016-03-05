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
    if(!collision_extent.norm()){
#if 0
        moi.translation=TV::Constant(1);
        moi.rotation=T_SPIN::Constant(1);
#else
        moi.translation=TV::Constant(linear_drag);
        moi.rotation=T_SPIN::Constant(rotational_drag);
#endif
    }
    else{
        T p=collision_extent.norm()/collision_radius+1;
        T E=sqrt(fabs(p*p-1))/p;
        T S=2*atanh(E)/E;
        
        // Perrin friction factors
        T half_length=collision_extent.norm()+collision_radius;
        T equivalent_linear_drag=6*M_PI*eta*std::pow(half_length*collision_radius*collision_radius,1.0f/3.0f);
        moi.translation=TV::Constant(2*std::pow(p,2.0f/3.0f)/S*equivalent_linear_drag);
        moi.rotation(2)=4.0f/3.0f*E*E/(2-S/(p*p))*rotational_drag;
        for(int spin=0;spin<2;spin++){
            moi.rotation(spin)=4.0f/3.0f*(1/(p*p)-p*p)/(2-S*(2-1/(p*p)))*rotational_drag;}}
    initialized=true;
}
///////////////////////////////////////////////////////////////////////
template<class TV> TV RIGID_STRUCTURE<TV>::
Displacement(const DATA<TV>& data,const RIGID_STRUCTURE<TV>& structure,TV& offset1,TV& offset2) const
{
    TV centroid1=frame.position;
    TV centroid2=centroid1+data.Minimum_Offset(frame.position,structure.frame.position);
    TV major_axis1=frame.orientation._transformVector(collision_extent);
    TV major_axis2=structure.frame.orientation._transformVector(structure.collision_extent);
    Matrix<T,2,1> weights;
    Matrix<TV,2,1> segment1;
    segment1[0]=centroid1-major_axis1;
    segment1[1]=centroid1+major_axis1;
    Matrix<TV,2,1> segment2;
    segment2[0]=centroid2-major_axis2;
    segment2[1]=centroid2+major_axis2;
    Segment_Segment_Displacement(segment1,segment2,weights);
    TV closest_point1=centroid1+(2*weights(0)-1)*major_axis1;
    TV closest_point2=centroid2+(2*weights(1)-1)*major_axis2;
    TV displacement=closest_point2-closest_point1;
    TV displacement_direction=displacement.normalized();
    offset1=closest_point1-centroid1;//+displacement_direction*collision_radius;
    offset2=closest_point2-centroid2;//+displacement_direction*structure->collision_radius;
    return displacement;
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RIGID_STRUCTURE)
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE,void)
{
    auto structure=std::make_shared<RIGID_STRUCTURE<TV>>();
    Parse_Vector(node["position"],structure->frame.position,TV());
    T scalar_collision_extent;
    Parse_Scalar(node["collision_extent"],scalar_collision_extent);
    structure->collision_extent=scalar_collision_extent*TV::UnitZ();
    Parse_Scalar(node["radius"],structure->radius,(T)0);
    Parse_Scalar(node["collision_radius"],structure->collision_radius,structure->radius);
    Parse_Scalar(node["kinematic"],structure->kinematic,false);
    Parse_Rotation(node["orientation"],structure->frame.orientation);
    structure->name=node["name"].asString();
    structure->Initialize_Inertia(simulation.data.globals["eta"]);
    auto rigid_structure_data=simulation.data.template Find_Or_Create<RIGID_STRUCTURE_DATA<TV>>();
    rigid_structure_data->structures.push_back(structure);
    LOG::cout<<"Name: "<<structure->name<<" index: "<<rigid_structure_data->structures.size()-1<<std::endl;
    return 0;
}
