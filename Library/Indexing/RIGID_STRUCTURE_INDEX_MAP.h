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
public:
    enum DEFINITIONS{STATIC_SIZE=TWIST<TV>::STATIC_SIZE,d=TV::RowsAtCompileTime};

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
        //unknown_map.template block<d,ROTATION<TV>::TwistSize>(0,d)=Cross_Product_Matrix(offset).transpose();
        return unknown_map;
    }

    static auto DP_DA(const TV& a,const TV& offset){
        auto norm_a=a.norm();
        auto s=sin(norm_a/2);
        auto c=cos(norm_a/2);
        auto dwda=-s/(2*norm_a)*a.transpose();
        auto q=s/2*a/norm_a;
        auto w=c;
        auto dqda=c/4*a*a.transpose()/sqr(norm_a)+s/2*(TV::Identity()/norm_a-a*a.transpose()/cube(norm_a));
        return -Cross_Product_Matrix(offset)*(2*dwda*q.transpose()+2*w*dqda);
    }

    static auto D2P_DA2(const TV& a,const TV& offset){
        auto norm_a=a.norm();
        auto aaT=a*a.transpose();
        auto s=sin(norm_a/2);
        auto c=cos(norm_a/2);
        auto dwda=-s/(2*norm_a)*a.transpose();
        auto q=s/2*a/norm_a;
        auto w=c;
        auto dqda=c/4*aaT/sqr(norm_a)+s/2*(TV::Identity()/norm_a-aaT/cube(norm_a));
        auto d2wda2=-dqda; // interestingly, this appears to be correct
        auto d2qda2=
        return -2*Cross_Product_Matrix(offset)*(d2wda2*q+2*dwda*dqda+w*d2qda2);
    }

};
}
#endif
