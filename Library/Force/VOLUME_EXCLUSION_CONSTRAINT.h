#ifndef __VOLUME_EXCLUSION_CONSTRAINT__
#define __VOLUME_EXCLUSION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class VOLUME_EXCLUSION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
public:
    struct CONSTRAINT{
        int s1;
        int s2;
    };
    std::vector<CONSTRAINT> constraints;

    VOLUME_EXCLUSION_CONSTRAINT(){}
    ~VOLUME_EXCLUSION_CONSTRAINT(){}

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    DEFINE_TYPE_NAME("VOLUME_EXCLUSION_CONSTRAINT")
};
}
#endif
