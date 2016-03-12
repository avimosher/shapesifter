#include <Analysis/ORIENTATION_QUALITY_PREDICATE.h>
#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar ORIENTATION_QUALITY_PREDICATE<TV>::
Sample_Point(const TV& point,const TV& offset,const DATA<TV>& data,RIGID_STRUCTURE<TV>& receptor,const RIGID_STRUCTURE<TV>& occluder,const T height,const T one_over_distance_limit_squared)
{
    receptor.frame.position=point+offset;
    TV o1,o2;
    if(occluder.Displacement(data,receptor,o1,o2).squaredNorm()>sqr(receptor.collision_radius+occluder.collision_radius)){
        return 1-(sqr(height)+offset.squaredNorm())*one_over_distance_limit_squared;}
    return 0;
}
//////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar ORIENTATION_QUALITY_PREDICATE<TV>::
Scalar(const SIMULATION<TV>& simulation)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto binder_structure=rigid_data->Structure(binder_name);
    const FRAME<TV>& binder_frame=binder_structure->frame;
    static const T two_thirds=2./3;
    static const T root_two_thirds=sqrt(two_thirds);
    static const T one_over_root_six=1/sqrt(6);
    static const T one_over_root_two=1/sqrt(2);
    TV binding_site=binder_frame*binder_bond_vector;
    TV o1,o2;
    T height=binding_site(2)-target_height;

    if(height<0 || height>distance_limit){return 0;}

    auto occluder=rigid_data->Structure(occluder_name);
    T distance_factor=(1-sqr(height/distance_limit))/2;
    if(occluder){
        TV central_point=binding_site;central_point(2)=occluder->frame.position(2);
        receptor.frame.position=central_point;
        TV direction=occluder->Displacement(simulation.data,receptor,o1,o2);
        T distance=direction.norm();direction.normalize();
        T contact_disk_radius=sqrt(sqr(distance_limit)-sqr(height));
        T one_over_distance_limit_squared=1/sqr(distance_limit);
        if(distance<contact_disk_radius+receptor.collision_radius+occluder->collision_radius){
            distance_factor=Sample_Point(central_point,TV(),simulation.data,receptor,*occluder,height,one_over_distance_limit_squared)/4;
            for(int s1=-1;s1<=1;s1+=2){
                distance_factor+=Sample_Point(central_point,s1*TV::Unit(0)*root_two_thirds*contact_disk_radius,simulation.data,receptor,*occluder,height,one_over_distance_limit_squared)/8;
                for(int s2=-1;s2<=1;s2+=2){
                    distance_factor+=Sample_Point(central_point,contact_disk_radius*(s1*TV::Unit(0)*one_over_root_six+s2*TV::Unit(2)*one_over_root_two),simulation.data,receptor,*occluder,height,one_over_distance_limit_squared)/8;}}}}

    T vertical_angle=Angle_Between(target_bond_orientation,binder_frame.orientation*binder_bond_vector);
    if(fabs(vertical_angle)>out_of_bond_angle_limit){return 0;}
    T vertical_angle_factor=1-sqr(vertical_angle/out_of_bond_angle_limit);
    T uniform_angle_factor=two_thirds/M_PI*around_bond_angle_limit;
    return distance_factor*vertical_angle_factor*uniform_angle_factor*simulation.dt;
}
//////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(ORIENTATION_QUALITY_PREDICATE)
DEFINE_AND_REGISTER_TEMPLATE_PARSER(ORIENTATION_QUALITY_PREDICATE,PREDICATE)
{
    auto predicate=std::make_shared<ORIENTATION_QUALITY_PREDICATE<TV>>();
    Parse_String(node["binder_name"],predicate->binder_name);
    Parse_String(node["occluder_name"],predicate->occluder_name);
    Parse_Scalar(node["distance_limit"],predicate->distance_limit);
    Parse_Scalar(node["around_bond_angle_limit"],predicate->around_bond_angle_limit);
    Parse_Scalar(node["out_of_bond_angle_limit"],predicate->out_of_bond_angle_limit);
    Parse_Vector(node["target_bond_orientation"],predicate->target_bond_orientation);
    Parse_Vector(node["binder_bond_vector"],predicate->binder_bond_vector);
    Parse_Scalar(node["target_height"],predicate->target_height);
    Parse_Scalar(node["receptor_radius"],predicate->receptor.collision_radius);
    return predicate;
}
