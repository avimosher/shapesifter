#ifndef __RELATIVE_POSITION_CONSTRAINT__
#define __RELATIVE_POSITION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
struct CONSTRAINT{
    typedef typename TV::Scalar T;
    int s1;
    TV v1;
    int s2;
    TV v2;
    T target_distance;
};

template<class TV>
class RELATIVE_POSITION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;

    std::vector<CONSTRAINT<TV>> constraints;
public:
    RELATIVE_POSITION_CONSTRAINT();
    ~RELATIVE_POSITION_CONSTRAINT();

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs);
    DEFINE_TYPE_NAME("RELATIVE_POSITION_CONSTRAINT")
};
}
#endif
