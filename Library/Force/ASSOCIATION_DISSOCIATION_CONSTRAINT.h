#ifndef __ASSOCIATION_DISSOCIATION_CONSTRAINT__
#define __ASSOCIATION_DISSOCIATION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class ASSOCIATION_DISSOCIATION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
public:
    using FORCE_TYPE<TV>::stored_forces;
    typedef int CONSTRAINT;
    std::unordered_map<CONSTRAINT,std::pair<int,T>> force_memory;

    struct INTERACTION{
        TV bond_vector;
        T bond_distance_threshold;
        T bond_orientation_threshold;

        T base_force_magnitude;
        T base_association_time;
        T base_dissociation_time;
        
        ROTATION<TV> relative_orientation;
        int s1;
        int s2;
    };

    std::vector<INTERACTION> interactions;
    std::vector<CONSTRAINT> constraints;

    ASSOCIATION_DISSOCIATION_CONSTRAINT(){}
    ~ASSOCIATION_DISSOCIATION_CONSTRAINT(){}

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("ASSOCIATION_DISSOCIATION_CONSTRAINT")
};
}
#endif
