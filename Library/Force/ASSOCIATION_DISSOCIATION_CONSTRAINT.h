#ifndef __ASSOCIATION_DISSOCIATION_CONSTRAINT__
#define __ASSOCIATION_DISSOCIATION_CONSTRAINT__

#include <Data/DATA.h>
#include <Data/ROTATION.h>
#include <Force/FORCE_TYPE.h>
#include <Utilities/HASHING.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>

namespace Mechanics{

template<class TV>
class PROXIMITY_SEARCH
{
    typedef typename TV::Scalar T;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime};
    const DATA<TV>& data;
    TV point;
    T proximity;
    std::vector<int> candidates;

    PROXIMITY_SEARCH(DATA<TV>& data_input,const TV& point_input,T proximity_input)
        :data(data_input),point(point_input),proximity(proximity_input)
    {}

    bool intersectVolume(const AlignedBox<T,d>& volume){
        TV half_edge_length=volume.sizes()/2;
        TV offset=data.Minimum_Offset(volume.center(),point);
        TV abs_offset=offset.cwiseAbs();
        TV difference=abs_offset-half_edge_length;
        return (difference.array()<=proximity).all();}

    bool intersectObject(int binder){
        candidates.push_back(binder);
        return false;}
};

template<class TV>
class STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT:public FORCE_REFERENCE<typename TV::Scalar>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    using FORCE_REFERENCE<T>::value;
    enum DEFINITIONS{t=T_SPIN::RowsAtCompileTime,d=TV::RowsAtCompileTime};
    typedef std::tuple<int,int,int> CONSTRAINT;
    std::vector<CONSTRAINT> constraints;


    STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT(){}
    virtual ~STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT(){}

    void setZero(){value.setZero();}
    virtual int Size(){return constraints.size()*(d+t);}
    DEFINE_TYPE_NAME("STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT");
};


template<class TV>
class ASSOCIATION_DISSOCIATION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::ORIENTATION T_ORIENTATION;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{t=T_SPIN::RowsAtCompileTime,d=TV::RowsAtCompileTime};
    typedef Matrix<T,d,t+d> LINEAR_CONSTRAINT_MATRIX;
    typedef Matrix<T,t,t+d> ANGULAR_CONSTRAINT_MATRIX;
    typedef Matrix<T,d+t,1> FORCE_VECTOR;
    using FORCE_TYPE<TV>::stored_forces;
    typedef std::tuple<int,int,int> CONSTRAINT;

    int call_count;
    std::unordered_map<CONSTRAINT,std::pair<int,FORCE_VECTOR>> force_memory;
    std::unordered_map<std::pair<int,int>,bool> partners;

    struct INTERACTION_TYPE{
        std::vector<std::vector<std::pair<int,bool>>> sites; // list of lists of body, active pairs
        std::vector<TV> site_offsets;
        T bond_distance_threshold;
        T bond_orientation_threshold;
        T base_association_time;
        T base_dissociation_time;
        ROTATION<TV> relative_orientation;
    };
    std::vector<INTERACTION_TYPE> interaction_types;
    std::vector<CONSTRAINT> constraints;

    ASSOCIATION_DISSOCIATION_CONSTRAINT()
        :call_count(1)
    {}
    ~ASSOCIATION_DISSOCIATION_CONSTRAINT(){}


    void Archive(cereal::BinaryOutputArchive& archive){
        archive(constraints,force_memory,partners,call_count);
        for(int i=0;i<interaction_types.size();i++){archive(interaction_types[i].sites);}
    }
    void Archive(cereal::BinaryInputArchive& archive){
        archive(constraints,force_memory,partners,call_count);
        for(int i=0;i<interaction_types.size();i++){archive(interaction_types[i].sites);}
    }

    std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const;
    void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information);
    void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_reference);
    void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_reference,int increment);
    ROTATION<TV> Find_Appropriate_Rotation(const ROTATION<TV>& rotation1,const ROTATION<TV>& rotation2);
    void Interaction_Candidates(DATA<TV>& data,const T dt,int type_index,int type1,int type2);
    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("ASSOCIATION_DISSOCIATION_CONSTRAINT")
};
}
#endif
