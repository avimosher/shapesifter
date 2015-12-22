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
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,1,t+d> CONSTRAINT_VECTOR;
    typedef Matrix<T,t+d,1> FORCE_VECTOR;
public:
    typedef std::tuple<int,int,int> CONSTRAINT; // body, axis, wall
 
    class STORED_WALL_CONSTRAINT:public FORCE_REFERENCE<T>
    {
    public:
        using FORCE_REFERENCE<T>::value;
        std::vector<CONSTRAINT> constraints;
        
        STORED_WALL_CONSTRAINT(){};
        virtual ~STORED_WALL_CONSTRAINT(){}
        
        void setZero(){value.setZero();}
        virtual int Size() const{return constraints.size();}
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

    virtual int DOF() const
    {return constraints.size();}

    std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const;
    void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment);
    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    DEFINE_TYPE_NAME("WALL_CONSTRAINT");
};
}
#endif
