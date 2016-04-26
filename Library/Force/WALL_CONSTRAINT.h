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
    enum MEMORY{MEMORY_COUNT,MEMORY_FORCE};
public:
    typedef std::array<int,4> INDICES; // body, substructure, axis, wall
 
    struct STORED_WALL_CONSTRAINT:public FORCE_REFERENCE<T>
    {
        using FORCE_REFERENCE<T>::value;
        std::vector<INDICES> constraints;
        
        STORED_WALL_CONSTRAINT(){};
        virtual ~STORED_WALL_CONSTRAINT(){}
        
        void setZero(){value.setZero();}
        virtual int Size() const{return constraints.size();}
        DEFINE_TYPE_NAME("STORED_WALL_CONSTRAINT");
    };

    struct CONSTRAINT
    {
        T_SPIN spin;
        T force;
        TV offset;
        FORCE_VECTOR force_direction;
        TV relative_position;

        CONSTRAINT(const T_SPIN& s,const T f,const TV& o,const FORCE_VECTOR& fd,const TV& rp)
            :spin(s),force(f),offset(o),force_direction(fd),relative_position(rp)
        {}
    };

    struct CONSTANT_FORCE
    {
        T_SPIN spin;
        TV offset;
        TV relative_position;
        T threshold;

        CONSTANT_FORCE(const T_SPIN& s,const TV& o,const TV& rp,const T t)
            :spin(s),offset(o),relative_position(rp),threshold(t)
        {}
    };

    RANGE<Matrix<bool,d,1>> walls;
    std::vector<CONSTRAINT> constraints;
    std::vector<INDICES> constraint_indices;
    std::vector<CONSTANT_FORCE> constant_forces;
    std::vector<INDICES> constant_force_indices;
    int call_count;
    std::unordered_map<INDICES,std::pair<int,T>> force_memory;
    std::unordered_map<INDICES,std::pair<int,T>> constant_force_memory;

    WALL_CONSTRAINT(){}
    ~WALL_CONSTRAINT(){}

    virtual int DOF() const
    {return constraints.size();}

    std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const;
    void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment);
    void Identify_Interactions_And_Compute_Errors(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    void Compute_Derivatives(DATA<TV>& data,FORCE<TV>& force,MATRIX_BUNDLE<TV>& system);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("WALL_CONSTRAINT");
};
}
#endif
