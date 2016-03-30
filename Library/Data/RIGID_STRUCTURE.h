#ifndef __RIGID_STRUCTURE__
#define __RIGID_STRUCTURE__

#include <Data/FRAME.h>
#include <Data/MOMENT.h>
#include <Data/TWIST.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class RIGID_STRUCTURE
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    bool initialized;
public:
    enum DEFINITIONS{PositionSize=FRAME<TV>::STATIC_SIZE,VelocitySize=TWIST<TV>::STATIC_SIZE};
    std::string name;
    FRAME<TV> frame;
    TWIST<TV> error;
    TWIST<TV> twist;
    MOMENT<TV> moi;
    T radius;
    T collision_radius;
    TV collision_extent; // require to be in the Z direction for now
    bool kinematic;
    std::shared_ptr<Matrix<T,4,1>> color;

    RIGID_STRUCTURE():initialized(false),kinematic(false){}
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


    DiagonalMatrix<T,VelocitySize> Inertia_Matrix()
    {
        assert(initialized);
        Matrix<T,VelocitySize,1> inertia;
        inertia<<moi.translation,
            moi.rotation;
        inertia*=-1;
        return inertia.asDiagonal();
    }

    AlignedBox<T,3> Bounding_Box(){
        TV minimum=frame*(-collision_extent);
        TV maximum=frame*collision_extent;
        return AlignedBox<T,3>(minimum.cwiseMin(maximum)-TV::Constant(radius),minimum.cwiseMax(maximum)+TV::Constant(radius));
    }

    template<class Archive>
    void serialize(Archive& archive) {archive(CEREAL_NVP(name),
            CEREAL_NVP(frame),
            CEREAL_NVP(error),
            CEREAL_NVP(moi),
            CEREAL_NVP(twist),
            CEREAL_NVP(radius),
            CEREAL_NVP(collision_radius),
            CEREAL_NVP(collision_extent),
            CEREAL_NVP(color),
            CEREAL_NVP(initialized));}

    TV Displacement(const DATA<TV>& data,const RIGID_STRUCTURE<TV>& structure,TV& offset1,TV& offset2) const;
    void Initialize_Inertia(const T eta);
    DEFINE_TYPE_NAME("RIGID_STRUCTURE")
};
}
#endif
