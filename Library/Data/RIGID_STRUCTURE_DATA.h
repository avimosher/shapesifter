#ifndef __RIGID_STRUCTURE_DATA__
#define __RIGID_STRUCTURE_DATA__

#include <Data/DATA_TYPE.h>
#include <Data/FRAME.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/TWIST.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/MATH.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <math.h>

namespace Mechanics{
template<class TV> class RIGID_STRUCTURE;

template<class TV>
class RIGID_STRUCTURE_DATA:public DATA_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{s=TWIST<TV>::STATIC_SIZE};
    std::vector<std::shared_ptr<RIGID_STRUCTURE<TV>>> structures;

    RIGID_STRUCTURE_DATA(){}
    ~RIGID_STRUCTURE_DATA(){}

    const Eigen::Rotation1D<T> Rotation_From_Angle_Axis(const Eigen::Matrix<T,0,1>& angle_axis) const
    {return Eigen::Rotation1D<T>();}

    const Eigen::Rotation2D<T> Rotation_From_Angle_Axis(const Eigen::Matrix<T,1,1>& angle_axis) const
    {return Rotation2D<T>(angle_axis.norm());}

    const Eigen::Quaternion<T> Rotation_From_Angle_Axis(const Eigen::Matrix<T,3,1>& rotation_vector) const {
        T spin_magnitude=rotation_vector.norm();
        Quaternion<T> q;
        q.w()=cos((T).5*spin_magnitude);
        q.vec()=(T).5*sinc((T).5*spin_magnitude)*rotation_vector;
        q.normalize();
        return q;
    }

    FRAME<TV> Updated_Frame(const DATA<TV>& data,const FRAME<TV>& frame,const TWIST<TV>& twist) const {
        FRAME<TV> updated_frame;updated_frame.position=data.Wrap(frame.position+twist.linear);
        updated_frame.orientation=Rotation_From_Angle_Axis(twist.angular)*frame.orientation;
        return updated_frame;
    }

    template<class Archive>
    void serialize(Archive& archive)
    {archive(structures);}

    int Structure_Index(const std::string& name) {
        for(int i=0;i<structures.size();i++){
            if(structures[i]->name==name){return i;}}
        return -1;
    }
    
    std::shared_ptr<RIGID_STRUCTURE<TV>> Structure(const std::string& name)
    {return structures[Structure_Index(name)];}

    int Size();
    virtual int Velocity_DOF() const;
    virtual int Position_DOF() const;
    virtual void Pack_Velocities(Block<Matrix<T,Dynamic,1>>& velocities);
    virtual void Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities);
    virtual void Pack_Positions(Block<Matrix<T,Dynamic,1>>& positions);
    virtual void Unpack_Positions(const Matrix<T,Dynamic,1>& positions);
    void Step(const DATA<TV>& data);
    virtual void Inertia(const T dt,std::vector<Triplet<T>>& force_terms,Matrix<T,Dynamic,1>& rhs);
    virtual void Viewer(osg::Node* node);

    DEFINE_TYPE_NAME("RIGID_STRUCTURE_DATA")
};
}
#endif
