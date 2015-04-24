#ifndef __ASSOCIATION_DISSOCIATION_CONSTRAINT__
#define __ASSOCIATION_DISSOCIATION_CONSTRAINT__

#include <Data/ROTATION.h>
#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class ASSOCIATION_DISSOCIATION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::ORIENTATION T_ORIENTATION;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{t=T_SPIN::RowsAtCompileTime,d=TV::RowsAtCompileTime};
    typedef Matrix<T,t,t> ROTATION_MATRIX;
    typedef Matrix<T,d,t+d> LINEAR_CONSTRAINT_MATRIX;
    typedef Matrix<T,t,t+d> ANGULAR_CONSTRAINT_MATRIX;
    typedef Matrix<T,t+d,1> FORCE_VECTOR;
    using FORCE_TYPE<TV>::stored_forces;
    typedef int CONSTRAINT;

    int call_count;
    std::unordered_map<CONSTRAINT,std::pair<int,FORCE_VECTOR>> force_memory;

    struct INTERACTION{
        T bond_distance_threshold;
        T bond_orientation_threshold;

        T base_association_time;
        T base_dissociation_time;
        
        ROTATION<TV> relative_orientation;
        TV v1;
        TV v2;
        int s1;
        int s2;
    };

    std::vector<INTERACTION> interactions;
    std::vector<CONSTRAINT> constraints;

    ASSOCIATION_DISSOCIATION_CONSTRAINT()
        :call_count(1)
    {}
    ~ASSOCIATION_DISSOCIATION_CONSTRAINT(){}

    void Unpack_Forces(const Matrix<T,Dynamic,1>& forces);
    void Increment_Forces(const Matrix<T,Dynamic,1>& forces);
    ROTATION<TV> Find_Appropriate_Rotation(const ROTATION<TV>& rotation1,const ROTATION<TV>& rotation2);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    DEFINE_TYPE_NAME("ASSOCIATION_DISSOCIATION_CONSTRAINT")
};
}
#endif
