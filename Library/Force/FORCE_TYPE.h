#ifndef __FORCE_TYPE__
#define __FORCE_TYPE__

#include <Utilities/TYPE_UTILITIES.h>
#include <Eigen/SparseCore>
#include <osg/Node>
#include <limits>

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

    int Force_DOF(){return stored_forces.rows();}
    virtual T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time){return std::numeric_limits<T>::max();}
    virtual void Pack_Forces(Block<Matrix<T,Dynamic,1>>& forces){forces=stored_forces;};
    virtual void Unpack_Forces(const Matrix<T,Dynamic,1>& forces){stored_forces=forces;};
    virtual void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)=0;
    template<class Archive> void serialize(Archive& archive){}
    virtual void Viewer(const DATA<TV>& data,osg::Node* node){};
    virtual std::string Name(){return "FORCE_TYPE";}
};
}
#endif
