#ifndef __ROTATION__
#define __ROTATION__

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace Mechanics{

template<class TV> struct ROTATION;

template<class T>
struct ROTATION<Eigen::Matrix<T,1,1>>
{
    enum DEFINITIONS{SizeAtCompileTime=0,TwistSize=0};
    typedef Eigen::Matrix<T,0,1> SPIN;
    typedef Eigen::Matrix<T,0,1> ORIENTATION;
};

template<class T>
struct ROTATION<Eigen::Matrix<T,2,1>>
{
    enum DEFINITIONS{SizeAtCompileTime=1,TwistSize=1};
    typedef Eigen::Matrix<T,1,1> SPIN;
    typedef Eigen::Rotation2D<T> ORIENTATION;
};

template<class T>
struct ROTATION<Eigen::Matrix<T,3,1>>
{
    enum DEFINITIONS{SizeAtCompileTime=4,TwistSize=3};
    typedef Eigen::Matrix<T,3,1> SPIN;
    typedef Eigen::Quaternion<T> ORIENTATION;
};


template<class T>
const Eigen::Matrix<T,4,1>& Get_Rotation_Coefficients(const Eigen::Quaternion<T>& rotation)
{return rotation.coeffs();}

template<class T>
Eigen::Matrix<T,1,1> Get_Rotation_Coefficients(const Eigen::Rotation2D<T>& rotation)
{
    Eigen::Matrix<T,1,1> matrix;matrix(0,0)=rotation.angle();return matrix;
}

template<class T>
const Eigen::Matrix<T,0,1>& Get_Rotation_Coefficients(const Eigen::Matrix<T,0,1>& rotation)
{return rotation;}

template<class T,class DERIVED>
void Set_Rotation_Coefficients(Eigen::Quaternion<T>& rotation,const Eigen::MatrixBase<DERIVED>& matrix)
{rotation.coeffs()=matrix;}

template<class T,class DERIVED>
void Set_Rotation_Coefficients(Eigen::Rotation2D<T>& rotation,const Eigen::MatrixBase<DERIVED>& matrix)
{rotation.angle()=matrix(0,0);}

template<class T,class DERIVED>
void Set_Rotation_Coefficients(Eigen::Matrix<T,0,1>& rotation,const Eigen::MatrixBase<DERIVED>& matrix)
{}



}

#endif

