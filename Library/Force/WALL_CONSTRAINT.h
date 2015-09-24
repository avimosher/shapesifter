#ifndef __WALL_CONSTRAINT__
#define __WALL_CONSTRAINT__

#include <Data/RANGE.h>
#include <Data/ROTATION.h>
#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class WALL_CONSTRAINT:public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef std::tuple<int,int,int> CONSTRAINT; // body, axis, wall
 
    class STORED_WALL_CONSTRAINT:public FORCE_REFERENCE<T>
    {
    public:
        using FORCE_REFERENCE<T>::value;
        std::vector<CONSTRAINT> constraints;
        
        STORED_WALL_CONSTRAINT(){};
        virtual ~STORED_WALL_CONSTRAINT(){}
        
        void setZero(){value.setZero();}
        virtual int Size(){return constraints.size();}
        DEFINE_TYPE_NAME("STORED_WALL_CONSTRAINT");
    };

   RANGE<Matrix<bool,d,1>> walls;

    std::vector<CONSTRAINT> constraints;
    std::vector<CONSTRAINT> constant_forces;
    int call_count;
    std::unordered_map<CONSTRAINT,std::pair<int,T>> force_memory;
    std::unordered_map<CONSTRAINT,std::pair<int,T>> constant_force_memory;

    WALL_CONSTRAINT()
    {}

    ~WALL_CONSTRAINT(){}

    std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const;
    void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    DEFINE_TYPE_NAME("WALL_CONSTRAINT");
};
}
#endif
