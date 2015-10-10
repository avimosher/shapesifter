#ifndef __FORCE_TYPE__
#define __FORCE_TYPE__

#include <Force/STORED_FORCE.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <limits>
#include <Eigen/SparseCore>
#include <osg/Node>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE;

template<class TV>
class FORCE_TYPE
{
    typedef typename TV::Scalar T;
public:
    Matrix<T,Dynamic,1> stored_forces;

    FORCE_TYPE(){}
    virtual ~FORCE_TYPE(){}

    virtual T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time){return std::numeric_limits<T>::max();}
    virtual std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const{return std::make_shared<FORCE_REFERENCE<T>>();}
    virtual void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information){force_information->value=stored_forces;};
    virtual void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information){stored_forces=force_information->value;};
    virtual void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment){stored_forces+=increment*force_information->value;};
    virtual void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)=0;
    virtual void Special(DATA<TV>& data,const T dt,const T time){};
    virtual void Archive(cereal::BinaryOutputArchive& archive){};
    virtual void Archive(cereal::BinaryInputArchive& archive){};
    virtual void Archive(cereal::JSONOutputArchive& archive){};
    virtual void Archive(cereal::JSONInputArchive& archive){};
    virtual void Viewer(const DATA<TV>& data,osg::Node* node){};
    virtual bool Equations_Changed() const{return false;}
    virtual std::string Name(){return "FORCE_TYPE";}
};
}
#endif
