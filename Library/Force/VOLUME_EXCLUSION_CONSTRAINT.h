#ifndef __VOLUME_EXCLUSION_CONSTRAINT__
#define __VOLUME_EXCLUSION_CONSTRAINT__

#include <Data/DATA.h>
#include <Force/FORCE_TYPE.h>
#include <Force/STORED_FORCE.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Utilities/CIRCULAR_STACK.h>

namespace Mechanics{
template<class TV> class RIGID_STRUCTURE;
typedef std::array<int,2> INDEX_PAIR;

template<class TV>
class BOX_PROXIMITY_SEARCH
{
    typedef typename TV::Scalar T;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime};
    const DATA<TV>& data;
    AlignedBox<T,d> box;
    std::vector<int> candidates;

    BOX_PROXIMITY_SEARCH(DATA<TV>& data_input,const AlignedBox<T,d>& input_box)
        :data(data_input),box(input_box)
    {}

    bool intersectVolume(const AlignedBox<T,d>& volume){
        TV offset=data.Minimum_Offset(volume.center(),box.center());
        TV abs_offset=offset.cwiseAbs();
        return (abs_offset.array()<=((volume.sizes()+box.sizes())/2).array()).all();}

    bool intersectObject(int structure){
        candidates.push_back(structure);
        return false;}
};

template<class T>
class STORED_VOLUME_EXCLUSION_CONSTRAINT:public FORCE_REFERENCE<T>
{
public:
    using FORCE_REFERENCE<T>::value;
    std::vector<INDEX_PAIR> constraints;

    STORED_VOLUME_EXCLUSION_CONSTRAINT(){};
    virtual ~STORED_VOLUME_EXCLUSION_CONSTRAINT(){}

    void setZero(){value.setZero();}
    virtual int Size() const{return constraints.size();}
    DEFINE_TYPE_NAME("STORED_VOLUME_EXCLUSION_CONSTRAINT");
};

template<class TV>
class VOLUME_EXCLUSION_CONSTRAINT:public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,1,t+d> CONSTRAINT_VECTOR;
    typedef Matrix<T,t+d,1> FORCE_VECTOR;
public:
    using FORCE_TYPE<TV>::stored_forces;using FORCE_TYPE<TV>::errors;


    enum CONSTRAINT_FIELDS{CONSTRAINT_FORCE_DIRECTIONS,CONSTRAINT_SPINS,CONSTRAINT_OFFSETS,CONSTRAINT_RELATIVE_POSITION};
    typedef std::tuple<std::array<FORCE_VECTOR,2>,std::array<T_SPIN,2>,std::array<TV,2>,TV> CONSTRAINT;
    /*struct CONSTRAINT{
        INDEX_PAIR indices;
        std::array<FORCE_VECTOR,2> force_directions;
        std::array<T_SPIN,2> spins;
        std::array<TV,2> offsets;
        TV relative_position;

        CONSTRAINT(const INDEX_PAIR& i_in,const std::array<FORCE_VECTOR,2>& fd_in,const std::array<T_SPIN,2>& s_in,const std::array<TV,2>& o_in,const TV& r_in)
            :indices(i_in),force_directions(fd_in),spins(s_in),offsets(o_in),relative_position(r_in)
        {}
        };*/

    enum CONSTANT_FORCE_FIELDS{CONSTANT_FORCE_SPINS,CONSTANT_FORCE_OFFSETS,CONSTANT_FORCE_RELATIVE_POSITION,CONSTANT_FORCE_THRESHOLD};
    typedef std::tuple<std::array<T_SPIN,2>,std::array<TV,2>,TV,T> CONSTANT_FORCE;

    /*struct CONSTANT_FORCE{
        INDEX_PAIR indices;
        std::array<T_SPIN,2> spins;
        std::array<TV,2> offsets;
        TV relative_position;
        T threshold;

        CONSTANT_FORCE(const INDEX_PAIR& i_in,const std::array<T_SPIN,2>& s_in,const std::array<TV,2>& o_in,const TV& r_in,const T t_in)
            :indices(i_in),spins(s_in),offsets(o_in),relative_position(r_in),threshold(t_in)
        {}
        };*/

    std::vector<INDEX_PAIR> constraint_indices;
    std::vector<INDEX_PAIR> constant_force_indices;
    std::vector<CONSTRAINT> constraints;
    std::vector<CONSTANT_FORCE> constant_forces;
    std::vector<T> rhs;
    int call_count;
    circular_stack<int> constraint_count;
    bool equations_changed;
    std::unordered_map<INDEX_PAIR,std::pair<int,T>> force_memory;
    std::unordered_map<INDEX_PAIR,std::pair<int,T>> constant_force_memory;

    VOLUME_EXCLUSION_CONSTRAINT()
        :call_count(0),constraint_count(2)
    {}
    ~VOLUME_EXCLUSION_CONSTRAINT(){}

    void Archive(cereal::BinaryOutputArchive& archive){archive(constraint_indices,constant_force_indices,force_memory,call_count,errors);}
    void Archive(cereal::BinaryInputArchive& archive){archive(constraint_indices,constant_force_indices,force_memory,call_count,errors);}

    void Archive(cereal::JSONOutputArchive& archive){archive(constraint_indices,constant_force_indices,force_memory,call_count,errors);}
    void Archive(cereal::JSONInputArchive& archive){archive(constraint_indices,constant_force_indices,force_memory,call_count,errors);}

    int DOF() const{return constraints.size();}
    std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const;
    void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment);
    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    void Identify_Interactions_And_Compute_Errors(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    void Compute_Derivatives(DATA<TV>& data,FORCE<TV>& force,MATRIX_BUNDLE<TV>& system);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    bool Equations_Changed() const{return equations_changed;};
    DEFINE_TYPE_NAME("VOLUME_EXCLUSION_CONSTRAINT")
};
}
#endif
