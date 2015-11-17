#ifndef DGP_VEC3_H
#define DGP_VEC3_H

#include <Eigen/Dense>
#include <Eigen/Geometry>


class Vec3:public Eigen::Matrix<double,3,1> {
public:
    Vec3( double X=0.0f, double Y=0.0f, double Z=0.0f ){
        Set(X,Y,Z);
    }

    Vec3(const Eigen::Matrix<double,3,1>& v){
        *this=v;
    }

    void Set( double X, double Y, double Z ){
        x()=X; y()=Y; z()=Z;
    }

    Vec3& operator=(const Eigen::Matrix<double,3,1>& v){
        Set(v[0],v[1],v[2]);
        return *this;
    }

    double& x(){return (*this)[0];}
    double& y(){return (*this)[1];}
    double& z(){return (*this)[2];}
};

#endif
