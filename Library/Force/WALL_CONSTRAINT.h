#ifndef __WALL_CONSTRAINT__
#define __WALL_CONSTRAINT__

namespace Mechanics{

template<class TV>
class WALL_CONSTRAINT:public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};

    enum WALL{PLUS_WALL,MINUS_WALL};
    typedef std::tuple<int,int,WALL> CONSTRAINT; // body, axis, wall
    RANGE<Matrix<bool,d>> walls;

    WALL_CONSTRAINT()
    {}

    ~WALL_CONSTRAINT(){}

    void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const;
    void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    DEFINE_TYPE_NAME("WALL_CONSTRAINT");
};
}
#endif
