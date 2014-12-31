#ifndef __RIGID_STRUCTURE_INDEX_MAP__
#define __RIGID_STRUCTURE_INDEX_MAP__

#include <Data/FRAME.h>

namespace Mechanics{
template<class TV>
class RIGID_STRUCTURE_INDEX_MAP
{
    typedef typename TV::Scalar T;
public:
    enum DEFINITIONS{STATIC_SIZE=FRAME<TV>::STATIC_SIZE};

    RIGID_STRUCTURE_INDEX_MAP();
    ~RIGID_STRUCTURE_INDEX_MAP();

    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Velocity_Map(const TV& offset){
        // map to TV
        // D rows, T+R columns.  D part is identity
        Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> map;
        for(int i=0;i<TV::RowsAtCompileTime;i++){map(i,i)=1;} // TODO: directly initialize to identity?
        Matrix<T,TV::RowsAtCompileTime,TV::RowsAtCompileTime> cross_product_matrix;
        cross_product_matrix<<0,-offset[2],offset[1],
            offset[2],0,-offset[0],
            -offset[1],offset[0],0;
        map.template block<3,3>(0,3)=cross_product_matrix;
    }

};
}
#endif
