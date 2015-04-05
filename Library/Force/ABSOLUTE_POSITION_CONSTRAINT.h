#ifndef __ABSOLUTE_POSITION_CONSTRAINT__
#define __ABSOLUTE_POSITION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

namespace Mechanics{

template<class TV>
class ABSOLUTE_POSITION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    using FORCE_TYPE<TV>::stored_forces;
    struct LINEAR_CONSTRAINT{
        int s;
        TV direction;
        T magnitude;

        template<class Archive>
        void serialize(Archive& archive)
        {archive(s,direction,magnitude);}
    };
    struct ANGULAR_CONSTRAINT{
        int s;
        ROTATION<TV> orientation;

        template<class Archive>
        void serialize(Archive& archive)
        {archive(s,orientation);}
    };
    std::vector<LINEAR_CONSTRAINT> linear_constraints;
    std::vector<ANGULAR_CONSTRAINT> angular_constraints;

    ABSOLUTE_POSITION_CONSTRAINT(){}
    ~ABSOLUTE_POSITION_CONSTRAINT(){}

    template<class Archive>
    void serialize(Archive& archive)
    {archive(linear_constraints,angular_constraints);}

    int Size()
    {return linear_constraints.size()+angular_constraints.size()*t;}

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("ABSOLUTE_POSITION_CONSTRAINT")
};
}
#endif
