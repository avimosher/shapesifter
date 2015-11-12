#ifndef __DEFINED_FORCE__
#define __DEFINED_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class DEFINED_FORCE:public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    struct FORCE_POKE{
        int s;
        TV v;
        TV f;
    };
    std::vector<FORCE_POKE> forces;

    DEFINED_FORCE(){}
    ~DEFINED_FORCE(){}

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    DEFINE_TYPE_NAME("DEFINED_FORCE")
};
}
#endif
