#ifndef __VOLUME_EXCLUSION_CONSTRAINT__
#define __VOLUME_EXCLUSION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class VOLUME_EXCLUSION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
public:
    using FORCE_TYPE<TV>::stored_forces;
    typedef std::pair<int,int> CONSTRAINT;
    std::vector<CONSTRAINT> constraints;
    int call_count;
    std::unordered_map<CONSTRAINT,std::pair<int,T>> force_memory;

    VOLUME_EXCLUSION_CONSTRAINT()
        :call_count(0)
    {}
    ~VOLUME_EXCLUSION_CONSTRAINT(){}

    void Unpack_Forces(const Matrix<T,Dynamic,1>& forces);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    DEFINE_TYPE_NAME("VOLUME_EXCLUSION_CONSTRAINT")
};
}
#endif
