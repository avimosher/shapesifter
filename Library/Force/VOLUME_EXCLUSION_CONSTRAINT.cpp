#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/VOLUME_EXCLUSION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <math.h>
using namespace Mechanics;
//#####################################################################
template<class TV> TV RADIAL_INEQUALITY_CONSTRAINT<TV>::
Segment_Closest_Points(const DATA<TV>& data,const FRAME<TV>& frame1,const FRAME<TV>& frame2,const RIGID_ATOM<TV>& rigid_atom1,const RIGID_ATOM<TV>& rigid_atom2,TV& world_space_offset1,TV& world_space_offset2)
{
    typedef typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT T_SEGMENT;
    TV centroid1=frame1.t;
    TV centroid2=centroid1+data.Minimum_Offset(frame1.t,frame2.t); // place so that closest distance algorithms will work

    TV major_axis1=frame1.r.Rotate(rigid_atom1.collision_extent);
    TV major_axis2=frame2.r.Rotate(rigid_atom2.collision_extent);
    T_SEGMENT segment1(centroid1-major_axis1,centroid1+major_axis1);
    T_SEGMENT segment2(centroid2-major_axis2,centroid2+major_axis2);

    VECTOR<T,2> weights;
    segment1.Shortest_Vector_Between_Segments(segment2,weights); // TODO: apparently not robust for parallel segments?  Check this.

    TV closest_point1=centroid1+(2*weights(1)-1)*major_axis1;
    TV closest_point2=centroid2+(2*weights(2)-1)*major_axis2;

    TV direction=closest_point2-closest_point1;

    world_space_offset1=closest_point1-centroid1+direction*rigid_atom1.collision_radius;
    world_space_offset2=closest_point2-centroid2+direction*rigid_atom2.collision_radius;
    return direction;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    T one_over_dt=1/dt;
    RANDOM<TV> random;
    
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.find("RIGID_STRUCTURE_DATA")->second);
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    RIGID_STRUCTURE_INDEX_MAP<TV> index_map;

    for(int s1=0;s1<rigid_data->structures.size();s1++){
        for(int s2=s1+1;s2<rigid_data->structures.size();s2++){
            

        }
    }

    constraint_rhs.resize(constraints.size(),1);
    for(int i=0;i<constraints.size();i++){
        const CONSTRAINT& constraint=constraints[i];
        int body_index1=constraint.s1;
        int body_index2=constraint.s2;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        FRAME<TV> frame1=rigid_structure1->frame;//rigid_data->Updated_Frame(data,rigid_structure1->frame,rigid_structure1->twist);
        FRAME<TV> frame2=rigid_structure2->frame;//rigid_data->Updated_Frame(data,rigid_structure2->frame,rigid_structure2->twist);
        TV offset1=frame1.orientation._transformVector(constraint.v1);
        TV offset2=frame2.orientation._transformVector(constraint.v2);

        TV direction=Segment_Closest_Points(data,frame1,frame2,rigid_structure1,rigid_structure2,offset1,offset2);
        T distance=direction.norm();
        T constraint_violation=distance-rigid_structure1->collision_radius-rigid_structure2->collision_radius;

        
        direction.normalize();
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index2,direction.transpose()*index_map.Velocity_Map(offset2)));
        terms.push_back(Triplet<CONSTRAINT_VECTOR>(i,body_index1,-direction.transpose()*index_map.Velocity_Map(offset1)));
        constraint_rhs(i,0)=constraint.target_distance-distance;
    }
    constraint_terms.resize(constraints.size(),RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE*rigid_data->structures.size());
    Flatten_Matrix(terms,constraint_terms);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(VOLUME_EXCLUSION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(VOLUME_EXCLUSION_CONSTRAINT,void)
{
    auto volume_exclusion_constraint=std::make_shared<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
    simulation.force.push_back(volume_exclusion_constraint);
    return 0;
}
