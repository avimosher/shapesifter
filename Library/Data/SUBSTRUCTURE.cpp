#include <Data/DATA.h>
#include <Data/SUBSTRUCTURE.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> TV SUBSTRUCTURE<TV>::
Displacement(const DATA<TV>& data,const FRAME<TV>& frame1,const FRAME<TV>& frame2,const SUBSTRUCTURE<TV>& substructure1,const SUBSTRUCTURE<TV>& substructure2,std::array<TV,2>& offsets)
{
    std::array<TV,2> centroids;
    centroids[0]=frame1.position+frame1.orientation*substructure1.offset;
    centroids[1]=centroids[0]+data.Minimum_Offset(centroids[0],frame2.position+frame2.orientation*substructure2.offset);
    std::array<TV,2> major_axes;
    major_axes[0]=frame1.orientation*substructure1.capsule_extent;
    major_axes[1]=frame2.orientation*substructure2.capsule_extent;
    Matrix<T,2,1> weights;
    std::array<Matrix<TV,2,1>,2> segments;
    for(int i=0;i<2;i++){
        segments[i][0]=centroids[i]-major_axes[i];
        segments[i][1]=centroids[i]+major_axes[i];}
    Segment_Segment_Displacement(segments,weights);
    std::array<TV,2> closest_points;
    for(int i=0;i<2;i++){
        closest_points[i]=centroids[i]+(2*weights(i)-1)*major_axes[i];
        offsets[i]=closest_points[i]-centroids[i];}
    return closest_points[1]-closest_points[0];
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SUBSTRUCTURE<TV>::
Parse(Json::Value& node)
{
    T scalar_capsule_extent;
    Parse_Scalar(node["capsule_extent"],scalar_capsule_extent);
    capsule_extent=scalar_capsule_extent*TV::UnitZ();
    Parse_Scalar(node["radius"],radius,(T)0);
    Parse_Vector(node["offset"],offset);
    if(!node["color"].isNull()){
        color.reset(new Matrix<T,4,1>);
        Parse_Vector(node["color"],*color);}
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(SUBSTRUCTURE)
GENERIC_TYPE_DEFINITION(SUBSTRUCTURE)
