#ifndef __RELATIVE_POSITION_CONSTRAINT__
#define __RELATIVE_POSITION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class RELATIVE_POSITION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    using FORCE_TYPE<TV>::stored_forces;
    struct CONSTRAINT{
        int s1;
        TV v1;
        int s2;
        TV v2;
        T target_distance;

        template<class Archive>
        void serialize(Archive& archive)
        {archive(s1,v1,s2,v2,target_distance);}
    };
    std::vector<CONSTRAINT> constraints;

    RELATIVE_POSITION_CONSTRAINT(){}
    ~RELATIVE_POSITION_CONSTRAINT(){}

    template<class Archive>
    void serialize(Archive& archive)
    {archive(constraints);}

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("RELATIVE_POSITION_CONSTRAINT")
};
}
#endif
