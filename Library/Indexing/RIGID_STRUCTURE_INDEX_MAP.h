#ifndef __RIGID_STRUCTURE_INDEX_MAP__
#define __RIGID_STRUCTURE_INDEX_MAP__

#include <Data/FRAME.h>
#include <Data/ROTATION.h>
#include <Utilities/EIGEN_HELPERS.h>

namespace Mechanics{
template<class TV>
class RIGID_STRUCTURE_INDEX_MAP
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{STATIC_SIZE=TWIST<TV>::STATIC_SIZE,d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};

    RIGID_STRUCTURE_INDEX_MAP(){}
    ~RIGID_STRUCTURE_INDEX_MAP(){}

    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Velocity_Map(const RIGID_STRUCTURE<TV>& structure,const TV& offset){
        return Velocity_Map(structure.frame.orientation._transformVector(offset));
    }

    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Velocity_Map(const TV& offset){
        // map to TV
        // D rows, T+R columns.  D part is identity
        Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> unknown_map;
        unknown_map.template block<d,d>(0,0).setIdentity();
        unknown_map.template block<d,ROTATION<TV>::TwistSize>(0,d)=Cross_Product_Matrix(offset).transpose();
        return unknown_map;
    }

    static Matrix<T,d,t+d> DD_DV(const RIGID_STRUCTURE<TV>& structure,const TV& object_offset,const TV& x1,const TV& x2,const TV& direction){
        T_SPIN a=structure.twist.angular;
        ROTATION<TV> orientation=structure.frame.orientation;
        orientation=ROTATION<TV>::From_Rotation_Vector(a).inverse()*orientation;
        TV offset=orientation._transformVector(object_offset);
        auto norm_a=a.norm();
        auto dw_da=-sinc(norm_a/2)/4*a.transpose();
        auto a_norma=(norm_a>1e-8)?(TV)(a/norm_a):TV::UnitX();
        Matrix<T,3,3> dq_da=cos(norm_a/2)/2*a_norma*a_norma.transpose()+sinc(norm_a/2)/2*(Matrix<T,t,t>::Identity()-a_norma*a_norma.transpose());
        //std::cout<<dq_da<<std::endl;
        auto w=cos(norm_a/2);
        auto q=sinc(norm_a/2)*a/2;
        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);
        //dx_da.template block<d,t>(0,d)=Cross_Product_Matrix(offset).transpose();

        std::cout<<"DX_DA: "<<dx_da<<std::endl;
        
        auto distance=std::max((T)1e-3,(x2-x1).norm());
        auto dd_da=dx_da/distance-direction/cube(distance)*direction.transpose()*dx_da;
        std::cout<<"DD_DA: "<<dd_da<<std::endl;
        return dd_da;
    }

    static Matrix<T,1,t+d> DC_DA(const RIGID_STRUCTURE<TV>& structure,const TV& object_offset,const TV& x1,const TV& x2,const TV& direction){
        T_SPIN a=structure.twist.angular;
        ROTATION<TV> orientation=structure.frame.orientation;
        orientation=ROTATION<TV>::From_Rotation_Vector(a).inverse()*orientation;
        TV offset=orientation._transformVector(object_offset);

        auto norm_a=a.norm();
        auto dw_da=-sinc(norm_a/2)/4*a.transpose();
        auto a_norma=(norm_a>1e-8)?(TV)(a/norm_a):TV::UnitX();
        Matrix<T,3,3> dq_da=cos(norm_a/2)/2*a_norma*a_norma.transpose()+sinc(norm_a/2)/2*(Matrix<T,t,t>::Identity()-a_norma*a_norma.transpose());
        //std::cout<<dq_da<<std::endl;
        auto w=cos(norm_a/2);
        auto q=sinc(norm_a/2)*a/2;
        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);
        //dx_da.template block<d,t>(0,d)=Cross_Product_Matrix(offset).transpose();

        //std::cout<<dx_da<<std::endl;
        
        auto distance=std::max((T)1e-3,(x2-x1).norm());
        auto dd_da=dx_da/distance-direction/cube(distance)*direction.transpose()*dx_da;
        //std::cout<<"DD_DA: "<<dd_da<<std::endl;
        TV normalized_direction=direction.normalized();
        //std::cout<<normalized_direction<<std::endl;
        auto final=direction.transpose()*dd_da+normalized_direction.transpose()*(dx_da);
        //auto final=normalized_direction.transpose()*(dx_da);
        //std::cout<<"DC_DA: "<<final<<std::endl;
        return final;
    }

    
};
}
#endif
