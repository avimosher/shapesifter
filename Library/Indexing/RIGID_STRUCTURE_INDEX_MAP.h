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
        // map to TV
        // D rows, T+R columns.  D part is identity
        Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> unknown_map;
        unknown_map.template block<d,d>(0,0).setIdentity();
        unknown_map.template block<d,ROTATION<TV>::TwistSize>(0,d)=Cross_Product_Matrix(structure.frame.rotation*offset).transpose();
        return unknown_map;
    }

};
}
#endif
