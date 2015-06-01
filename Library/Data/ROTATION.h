#ifndef __ROTATION__
#define __ROTATION__

#include <Utilities/EIGEN_HELPERS.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// TODO: move to math Utilities header
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

namespace Mechanics{

template<class TV> class ROTATION;

template<class T>
class ROTATION<Eigen::Matrix<T,1,1>>:public Eigen::Rotation1D<T>
{
    typedef Eigen::Matrix<T,0,1> TV;
  public:
    enum DEFINITIONS{SizeAtCompileTime=0,TwistSize=0};
    typedef Eigen::Matrix<T,0,1> SPIN;
    typedef Eigen::Rotation1D<T> ORIENTATION;

    ROTATION(){}
    ROTATION(const Eigen::Rotation1D<T>& r)
        :Eigen::Rotation1D<T>(r)
    {}

    TV Axis() const {return TV();}
    T Angle() const {return T();}
    ROTATION Scale_Angle(T scale){return *this;}
    int Sign() const{return 1;}
    TV Vec() const{return TV();}
    T W() const{return 0;}
    static ROTATION From_Rotation_Vector(const TV&){return ROTATION(Eigen::Rotation1D<T>());}
};

template<class T>
class ROTATION<Eigen::Matrix<T,2,1>>:public Eigen::Rotation2D<T>
{
    typedef Eigen::Matrix<T,1,1> TV;
  public:
    enum DEFINITIONS{SizeAtCompileTime=1,TwistSize=1};
    typedef Eigen::Matrix<T,1,1> SPIN;
    typedef Eigen::Rotation2D<T> ORIENTATION;

    ROTATION(){}
    ROTATION(const Eigen::Rotation2D<T>& r)
        :Eigen::Rotation2D<T>(r)
    {}

    TV Axis() const {return TV();}
    T Angle() const {return (*this).angle();}
    ROTATION Scale_Angle(T scale){
        (*this).angle()*=scale;
        return *this;
    }
    int Sign() const{return 1;}
    TV Vec() const{return TV();}
    T W() const{return 0;}
    static ROTATION From_Rotation_Vector(const TV&){return ROTATION(0);}
};

template<class T>
class ROTATION<Eigen::Matrix<T,3,1>>:public Eigen::Quaternion<T>
{
    typedef Eigen::Matrix<T,3,1> TV;
  public:
    enum DEFINITIONS{SizeAtCompileTime=4,TwistSize=3};
    typedef Eigen::Matrix<T,3,1> SPIN;
    typedef ROTATION ORIENTATION;

    ROTATION()
    {
        Eigen::Quaternion<T>::setIdentity();
    }
    ROTATION(const Eigen::Quaternion<T>& q)
        :Eigen::Quaternion<T>(q)
    {}

    ROTATION operator*(const ROTATION& rhs) const
    {return ((Eigen::Quaternion<T>)(*this))*((Eigen::Quaternion<T>)rhs);}

    TV operator*(const TV& rhs) const
    {return this->_transformVector(rhs);}

    TV Axis() const {return Eigen::AngleAxis<T>(*this).axis();}
    T Angle() const {return Eigen::AngleAxis<T>(*this).angle();}

    ROTATION<TV> Scale_Angle(T scale){
        Eigen::AngleAxis<T> angle_axis(*this);angle_axis.angle()*=scale;
        (*this)=Eigen::Quaternion<T>(angle_axis);
        return *this;
    }

    int Sign() const{return sgn((*this).w());}
    TV Vec() const{return (*this).vec();}
    T W() const{return (*this).w();}
    TV Rotation_Vector() const{return Angle()*Axis();}
    static ROTATION From_Rotation_Vector(const TV& rotation_vector){
        T norm=rotation_vector.norm();
        return Eigen::Quaternion<T>(Eigen::AngleAxis<T>(norm,norm?rotation_vector.normalized():TV::UnitX()));
    }

    static ROTATION From_Rotated_Vector(const TV& initial_vector,const TV& final_vector){
        return Eigen::Quaternion<T>::FromTwoVectors(initial_vector,final_vector);
    }

    static ROTATION From_Rotated_Vector_Around_Axis(const TV& initial_vector,const TV& final_vector,const TV& axis_direction){
        TV axis_normalized=axis_direction.normalized();
        T v0_projection=initial_vector.dot(axis_normalized);
        TV v0_base=initial_vector-v0_projection*axis_normalized;
        T v1_projection=final_vector.dot(axis_normalized);
        TV v1_base=final_vector-v1_projection*axis_normalized;
        TV v0=v0_base.normalized();
        TV v1=v1_base.normalized();
        
        TV axis=v0.cross(v1);
        Eigen::Quaternion<T> q;
        q.w()=1+v0.dot(v1);
        q.vec()=axis;
        q.normalize();
        return q;
    }
};

template<class T>
const Eigen::Matrix<T,4,1>& Get_Rotation_Coefficients(const Eigen::Quaternion<T>& rotation)
{return rotation.coeffs();}

template<class T>
Eigen::Matrix<T,1,1> Get_Rotation_Coefficients(const Eigen::Rotation2D<T>& rotation)
{Eigen::Matrix<T,1,1> matrix;matrix(0,0)=rotation.angle();return matrix;}

template<class T>
Eigen::Matrix<T,0,1> Get_Rotation_Coefficients(const Eigen::Rotation1D<T>& rotation)
{return Eigen::Matrix<T,0,1>();}

template<class T,class DERIVED>
void Set_Rotation_Coefficients(Eigen::Quaternion<T>& rotation,const Eigen::MatrixBase<DERIVED>& matrix)
{rotation.coeffs()=matrix;}

template<class T,class DERIVED>
void Set_Rotation_Coefficients(Eigen::Rotation2D<T>& rotation,const Eigen::MatrixBase<DERIVED>& matrix)
{rotation.angle()=matrix(0,0);}

template<class T,class DERIVED>
void Set_Rotation_Coefficients(Eigen::Rotation1D<T>& rotation,const Eigen::MatrixBase<DERIVED>& matrix)
{}

}

#endif

