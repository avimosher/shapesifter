#ifndef __ROTATION__
#define __ROTATION__

#include <Utilities/EIGEN_HELPERS.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace Mechanics{

template<class TV> class ROTATION;

template<class T>
class ROTATION<Eigen::Matrix<T,1,1>>:public Eigen::Rotation1D<T>
{
  public:
    enum DEFINITIONS{SizeAtCompileTime=0,TwistSize=0};
    typedef Eigen::Matrix<T,0,1> SPIN;
    typedef Eigen::Rotation1D<T> ORIENTATION;
};

template<class T>
class ROTATION<Eigen::Matrix<T,2,1>>:public Eigen::Rotation2D<T>
{
  public:
    enum DEFINITIONS{SizeAtCompileTime=1,TwistSize=1};
    typedef Eigen::Matrix<T,1,1> SPIN;
    typedef Eigen::Rotation2D<T> ORIENTATION;
};

template<class T>
class ROTATION<Eigen::Matrix<T,3,1>>:public Eigen::Quaternion<T>
{
    typedef Eigen::Matrix<T,3,1> TV;
  public:
    enum DEFINITIONS{SizeAtCompileTime=4,TwistSize=3};
    typedef Eigen::Matrix<T,3,1> SPIN;
    typedef Eigen::Quaternion<T> ORIENTATION;

    TV operator*(const TV& vector){
        return _transformVector(vector);
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

