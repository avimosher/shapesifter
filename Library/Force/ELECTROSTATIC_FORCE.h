#ifndef __ELECTROSTATIC_FORCE__
#define __ELECTROSTATIC_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class ELECTROSTATIC_FORCE : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};

    typedef std::tuple<int,T,TV> CHARGE; // body, charge, offset
    std::vector<CHARGE> charges;

    ELECTROSTATIC_FORCE(){}
    ~ELECTROSTATIC_FORCE(){}

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);

    void Parse_Fields();
    DEFINE_TYPE_NAME("ELECTROSTATIC_FORCE");
};
}
#endif
