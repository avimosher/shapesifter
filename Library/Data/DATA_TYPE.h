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
    virtual int Velocity_DOF() const {return 0;}
    virtual int Position_DOF() const {return 0;}
    virtual void Pack_Velocities(Block<Matrix<T,Dynamic,1>>& velocities){};
    virtual void Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities){};
    virtual void Pack_Positions(Block<Matrix<T,Dynamic,1>>& positions){};
    virtual void Unpack_Positions(const Matrix<T,Dynamic,1>& positions){};
    virtual void Step(const DATA<TV>& data){};
    template<class Archive> void serialize(Archive& archive){}
    virtual T Print(){return T();}
    virtual void Viewer(osg::Node* node){};
    virtual std::string Name(){return "DATA_TYPE";}
    virtual void Inertia(const T dt,std::vector<Triplet<T>>& force_terms,Matrix<T,Dynamic,1>& rhs){};
};

}
#endif
