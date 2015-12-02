#ifndef __DATA_TYPE__
#define __DATA_TYPE__

#include <Utilities/TYPE_UTILITIES.h>
#include <osg/Node>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class DATA_TYPE
{
    typedef typename TV::Scalar T;

public:
    DATA_TYPE(){}
    ~DATA_TYPE(){}

    virtual int Size()=0;
    virtual int DOF() const{return Velocity_DOF();}
    virtual int Velocity_DOF() const {return 0;}
    virtual int Position_DOF() const {return 0;}
    virtual void Identify_DOF(int index) const{LOG::cout<<Name()<<" DOF "<<index<<std::endl;}
    virtual void Pack_Velocities(Block<Matrix<T,Dynamic,1>>& velocities){};
    virtual void Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities){};
    virtual void Pack_Positions(Block<Matrix<T,Dynamic,1>>& positions){};
    virtual void Unpack_Positions(const Matrix<T,Dynamic,1>& positions){};
    virtual void Store_Errors(const Matrix<T,Dynamic,1>& errors){};
    virtual void Step(const DATA<TV>& data){};
    virtual void Eliminate_Rotation(const DATA<TV>& data){};
    template<class Archive> void serialize(Archive& archive){}
    virtual T Print(){return T();}
    virtual void Viewer(osg::Node* node){};
    virtual std::string Name() const{return "DATA_TYPE";}
    virtual void Inertia(const T dt,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& inverse_inertia,Matrix<T,Dynamic,1>& rhs){};
    virtual void Kinematic_Projection(SparseMatrix<T>& kinematic_projection){};
};

}
#endif
