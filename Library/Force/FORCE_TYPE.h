#ifndef __FORCE_TYPE__
#define __FORCE_TYPE__

#include <Force/STORED_FORCE.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <limits>
#include <Eigen/SparseCore>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE;
template<class TV> class MATRIX_BUNDLE;

template<class TV>
class FORCE_TYPE
{
    typedef typename TV::Scalar T;
public:
    Matrix<T,Dynamic,1> stored_forces;
    Matrix<T,Dynamic,1> errors;

    FORCE_TYPE(){}
    virtual ~FORCE_TYPE(){}

    virtual T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time){return std::numeric_limits<T>::max();}
    virtual std::shared_ptr<FORCE_REFERENCE<T>> Create_Stored_Force() const{return std::make_shared<FORCE_REFERENCE<T>>();}
    virtual void Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information){force_information->value=stored_forces;};
    virtual void Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information){stored_forces=force_information->value;};
    virtual void Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment){stored_forces+=increment*force_information->value;}
    virtual void Store_Errors(std::shared_ptr<FORCE_REFERENCE<T>> force_information){errors=force_information->value;};
    
    virtual void Identify_Interactions_And_Compute_Errors(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,MATRIX_BUNDLE<TV>& system,bool stochastic)
    {Linearize(data,force,dt,time,system,stochastic);}
    virtual void Compute_Derivatives(DATA<TV>& data,FORCE<TV>& force,MATRIX_BUNDLE<TV>& system){};
    virtual void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time,MATRIX_BUNDLE<TV>& system,bool stochastic){};
    virtual int DOF() const{return stored_forces.size();}
    virtual int Rows() const{return DOF();}
    virtual int Columns() const{return DOF();}
    virtual void Identify_DOF(int index) const{LOG::cout<<Name()<<" DOF "<<index<<std::endl;}
    virtual void Archive(cereal::BinaryOutputArchive& archive){};
    virtual void Archive(cereal::BinaryInputArchive& archive){};
    virtual void Archive(cereal::JSONOutputArchive& archive){};
    virtual void Archive(cereal::JSONInputArchive& archive){};
    virtual bool Equations_Changed() const{return false;}
    virtual std::string Name() const{return "FORCE_TYPE";}
};
}
#endif
