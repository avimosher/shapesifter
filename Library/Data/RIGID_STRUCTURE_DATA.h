//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class RIGID_STRUCTURE_DATA
///////////////////////////////////////////////////////////////////////
#ifndef __RIGID_STRUCTURE_DATA__
#define __RIGID_STRUCTURE_DATA__

#include <Data/DATA_TYPE.h>
#include <Data/FRAME.h>
#include <Data/RIGID_STRUCTURE.h>
#include <Data/TWIST.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

namespace Mechanics{
template<class TV> class RIGID_STRUCTURE;

template<class TV>
class RIGID_STRUCTURE_DATA:public DATA_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    std::vector<std::shared_ptr<RIGID_STRUCTURE<TV>>> structures;

    RIGID_STRUCTURE_DATA();
    ~RIGID_STRUCTURE_DATA(){}

    const Eigen::Rotation1D<T> Rotation_From_Angle_Axis(const Eigen::Matrix<T,0,1>& angle_axis) const {
        return Eigen::Rotation1D<T>();
    }

    const Eigen::Rotation2D<T> Rotation_From_Angle_Axis(const Eigen::Matrix<T,1,1>& angle_axis) const {
        return Rotation2D<T>(angle_axis.norm());
    }

    const Eigen::AngleAxis<T> Rotation_From_Angle_Axis(const Eigen::Matrix<T,3,1>& angle_axis) const {
        Eigen::Matrix<T,3,1> spin_axis=angle_axis.normalized();
        T spin_magnitude=angle_axis.norm();
        return Eigen::AngleAxis<T>(spin_magnitude,spin_axis);
    }

    FRAME<TV> Updated_Frame(const DATA<TV>& data,const FRAME<TV>& frame,const TWIST<TV>& twist) const {
        FRAME<TV> updated_frame;updated_frame.position=data.Wrap(frame.position+twist.linear);
        updated_frame.orientation=Rotation_From_Angle_Axis(twist.angular)*frame.orientation;
        return updated_frame;
    }

    template<class Archive>
    void serialize(Archive& archive) {
        //archive(framing);
        archive(structures);
    }

    int Size();
    Matrix<T,Dynamic,1> Variables();
    void Step(const Matrix<T,Dynamic,1>& variables);
    virtual void Viewer(osg::Node* node);

    DEFINE_TYPE_NAME("RIGID_STRUCTURE_DATA")
};
}
#endif
