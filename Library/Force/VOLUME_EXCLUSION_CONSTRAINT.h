#ifndef __VOLUME_EXCLUSION_CONSTRAINT__
#define __VOLUME_EXCLUSION_CONSTRAINT__

#include <Data/DATA.h>
#include <Force/FORCE_TYPE.h>
#include <Force/STORED_FORCE.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Utilities/CEREAL_HELPERS.h>

namespace Mechanics{
template<class TV> class RIGID_STRUCTURE;
typedef std::pair<int,int> CONSTRAINT;

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
    std::vector<CONSTRAINT> constraints;

    STORED_VOLUME_EXCLUSION_CONSTRAINT(){};
    virtual ~STORED_VOLUME_EXCLUSION_CONSTRAINT(){}

    void setZero(){value.setZero();}
    virtual int Size(){return constraints.size();}
    DEFINE_TYPE_NAME("STORED_VOLUME_EXCLUSION_CONSTRAINT");
};

template<class TV>
class VOLUME_EXCLUSION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    using FORCE_TYPE<TV>::stored_forces;
    std::vector<CONSTRAINT> constraints;
    std::vector<CONSTRAINT> constant_forces;
    int call_count;
    std::unordered_map<CONSTRAINT,std::pair<int,T>> force_memory;
    std::unordered_map<CONSTRAINT,std::tuple<int,T,T,T,T>> constant_force_memory;

    VOLUME_EXCLUSION_CONSTRAINT()
        :call_count(0)
    {}
    ~VOLUME_EXCLUSION_CONSTRAINT(){}

    void Archive(cereal::BinaryOutputArchive& archive){archive(constraints,constant_forces,force_memory,call_count);}
    void Archive(cereal::BinaryInputArchive& archive){archive(constraints,constant_forces,force_memory,call_count);}

    void Archive(cereal::JSONOutputArchive& archive){archive(constraints,constant_forces,force_memory,call_count);}
    void Archive(cereal::JSONInputArchive& archive){archive(constraints,constant_forces,force_memory,call_count);}

    int Force_DOF(){return constraints.size();}
    std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const;
    void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("VOLUME_EXCLUSION_CONSTRAINT")
};
}
#endif
