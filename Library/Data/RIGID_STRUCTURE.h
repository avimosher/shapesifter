#ifndef __RIGID_STRUCTURE__
#define __RIGID_STRUCTURE__

#include <Data/FRAME.h>
#include <Data/TWIST.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <cereal/archives/binary.hpp>

namespace Mechanics{

template<class TV>
class RIGID_STRUCTURE
{
    typedef typename TV::Scalar T;
public:
    enum DEFINITIONS{STATIC_SIZE=FRAME<TV>::STATIC_SIZE};
    std::string name;
    FRAME<TV> frame;
    TWIST<TV> twist;
    T radius;
    T collision_radius;
    T collision_extent; // defined in the Z direction
#if 0
    MOMENT<TV> moi;
#endif

    RIGID_STRUCTURE(){}
    ~RIGID_STRUCTURE(){}

    static Matrix<T,1,1> Segment_Segment_Displacement(const Matrix<Matrix<T,1,1>,2,1>& s1,const Matrix<Matrix<T,1,1>,2,1>& s2,Matrix<T,2,1>& weights){return Matrix<T,1,1>();}

    static Matrix<T,2,1> Segment_Segment_Displacement(const Matrix<Matrix<T,2,1>,2,1>& s1,const Matrix<Matrix<T,2,1>,2,1>& s2,Matrix<T,2,1>& weights){return Matrix<T,2,1>();}

    // From PhysBAM
    static Matrix<T,3,1> Segment_Segment_Displacement(const Matrix<Matrix<T,3,1>,2,1>& s1,const Matrix<Matrix<T,3,1>,2,1>& s2,Matrix<T,2,1>& weights)
    {
        Matrix<T,3,1> u=s1[1]-s1[0],v=s2[1]-s2[0],w=s2[0]-s1[0];
        T u_magnitude_squared=u.squaredNorm(),
            v_magnitude_squared=v.squaredNorm(),
            u_dot_u=u_magnitude_squared,
            v_dot_v=v_magnitude_squared,
            u_dot_v=u.dot(v),
            u_dot_w=u.dot(w),
            v_dot_w=v.dot(w);
        T v_dot_v_pseudoinverse=(v_dot_v?1/v_dot_v:0);
        T denominator=u_dot_u*v_dot_v-u_dot_v*u_dot_v,rhs1=v_dot_v*u_dot_w-u_dot_v*v_dot_w,rhs2=u_dot_v*u_dot_w-u_dot_u*v_dot_w;
        bool check_boundary=false;
        if(rhs1<=0 || denominator<=rhs1){check_boundary=true;}
        else{weights[0]=rhs1/denominator;}
        if(rhs2<=0 || denominator<=rhs2){check_boundary=true;}
        else{weights[1]=rhs2/denominator;}
        if(check_boundary){
            T v_plus_w_dot_u=u_dot_v+u_dot_w,u_minus_w_dot_v=u_dot_v-v_dot_w,distance_squared_minus_w_dot_w;
            weights[0]=0; // check weights[0]=0 side

            if(v_dot_w>=0){distance_squared_minus_w_dot_w=0;weights[1]=0;}
            else if(v_dot_v<=-v_dot_w){distance_squared_minus_w_dot_w=v_dot_v+2*v_dot_w;weights[1]=1;}
            else{weights[1]=-v_dot_w/v_dot_v;distance_squared_minus_w_dot_w=weights[1]*v_dot_w;}
            // check weights.x=1 side
            if(u_minus_w_dot_v<=0){T new_distance_squared=u_dot_u-2*u_dot_w;
                if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights<<1,0;}}
            else if(v_dot_v<=u_minus_w_dot_v){T new_distance_squared=v_dot_v+2*(v_dot_w-v_plus_w_dot_u)+u_dot_u;
                if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights<<1,1;}}
            else{T weights_y_temp=u_minus_w_dot_v*v_dot_v_pseudoinverse,new_distance_squared=u_dot_u-2*u_dot_w-weights_y_temp*u_minus_w_dot_v;
                if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights<<1,weights_y_temp;}}
            // check weights.y=0 side ignoring corners (already handled above)
            if(u_dot_w>0 && u_dot_u>u_dot_w){T weights_x_temp=u_dot_w/u_dot_u,new_distance_squared=-weights_x_temp*u_dot_w;
                if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights<<weights_x_temp,0;}}
            // check weights.y=1 side ignoring corners (already handled above)
            if(v_plus_w_dot_u>0 && u_dot_u>v_plus_w_dot_u){T weights_x_temp=v_plus_w_dot_u/u_dot_u,new_distance_squared=v_dot_v+2*v_dot_w-weights_x_temp*v_plus_w_dot_u;
                if(new_distance_squared<distance_squared_minus_w_dot_w){distance_squared_minus_w_dot_w=new_distance_squared;weights<<weights_x_temp,1;}}}
        return weights[0]*u-w-weights[1]*v;
    }

    TV Displacement(const DATA<TV>& data,const std::shared_ptr<RIGID_STRUCTURE<TV>> structure,TV& offset1,TV& offset2) const {
        TV centroid1=frame.position;
        TV centroid2=centroid1+data.Minimum_Offset(frame.position,structure->frame.position);
        TV major_axis1=frame.orientation._transformVector(collision_extent*TV::UnitZ());
        TV major_axis2=structure->frame.orientation._transformVector(structure->collision_extent*TV::UnitZ());
        Matrix<T,2,1> weights;
        Matrix<TV,2,1> segment1;
        segment1[0]=centroid1-major_axis1;
        segment1[1]=centroid1+major_axis1;
        Matrix<TV,2,1> segment2;
        segment2[0]=centroid2-major_axis2;
        segment2[1]=centroid2+major_axis2;
        Segment_Segment_Displacement(segment1,segment2,weights);
        TV closest_point1=centroid1+(2*weights(1)-1)*major_axis1;
        TV closest_point2=centroid2+(2*weights(1)-1)*major_axis2;
        TV displacement=closest_point2-closest_point1;
        offset1=closest_point1-centroid1+displacement*collision_radius;
        offset2=closest_point2-centroid2+displacement*structure->collision_radius;
        return displacement;
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(name,frame,twist,radius,collision_radius,collision_extent);
    }

    DEFINE_TYPE_NAME("RIGID_STRUCTURE")
};
}
#endif
