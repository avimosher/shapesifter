#ifndef __RELATIVE_POSITION_CONSTRAINT__
#define __RELATIVE_POSITION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV>
class RELATIVE_POSITION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;

public:
    T temperature;

    RELATIVE_POSITION_CONSTRAINT();
    ~RELATIVE_POSITION_CONSTRAINT();

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side);    
    DEFINE_TYPE_NAME("RELATIVE_POSITION_CONSTRAINT")
};
}
#endif
