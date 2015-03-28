#ifndef __RIGID_STRUCTURE_INDEX_MAP__
#define __RIGID_STRUCTURE_INDEX_MAP__

#include <Data/FRAME.h>
#include <Data/ROTATION.h>
#include <Utilities/EIGEN_HELPERS.h>

namespace Mechanics{
template<class TV>
class RIGID_STRUCTURE_INDEX_MAP
{
    typedef typename TV::Scalar T;
    typedef Matrix<T,1,TV::RowsAtCompileTime> TV_T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{STATIC_SIZE=TWIST<TV>::STATIC_SIZE,d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};

    RIGID_STRUCTURE_INDEX_MAP(){}
    ~RIGID_STRUCTURE_INDEX_MAP(){}

    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Velocity_Map(const RIGID_STRUCTURE<TV>& structure,const TV& offset){
        return Velocity_Map(structure.frame.orientation*offset);
    }

    static Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> Velocity_Map(const TV& offset){
        // map to TV
        // D rows, T+R columns.  D part is identity
        Matrix<T,TV::RowsAtCompileTime,STATIC_SIZE> unknown_map;
        unknown_map.template block<d,d>(0,0).setIdentity();
        unknown_map.template block<d,ROTATION<TV>::TwistSize>(0,d)=Cross_Product_Matrix(offset).transpose();
        return unknown_map;
    }

    static Matrix<T,t+d,t+d> DF_DA(const RIGID_STRUCTURE<TV>& structure,const TV& object_offset,const TV& x1,const TV& x2,const TV& direction){
        T_SPIN a=structure.twist.angular;
        ROTATION<TV> orientation=structure.frame.orientation;
        TV current_offset=orientation*object_offset;
        orientation=ROTATION<TV>::From_Rotation_Vector(a).inverse()*orientation;
        TV offset=orientation*object_offset;
        T norm_a=a.norm();
        TV_T dw_da=-sinc(norm_a/2)/4*a.transpose();
        TV a_norma=(norm_a>1e-8)?(TV)(a/norm_a):TV::UnitX();
        Matrix<T,3,3> dq_da=cos(norm_a/2)/2*a_norma*a_norma.transpose()+sinc(norm_a/2)/2*(Matrix<T,t,t>::Identity()-a_norma*a_norma.transpose());
        T w=cos(norm_a/2);
        TV q=sinc(norm_a/2)*a/2;
        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);
        //dx_da.template block<d,t>(0,d)=Cross_Product_Matrix(offset).transpose();

        T distance=std::max((T)1e-3,(x2-x1).norm());
        Matrix<T,d,t+d> dd_da=dx_da/distance-direction/cube(distance)*direction.transpose()*dx_da;

        TV d_n=direction.normalized();
        Matrix<T,d,t+d> dr_da;
        dr_da.template block<d,d>(0,0).setZero();
        dr_da.template block<d,t>(0,d)=-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);

        Matrix<T,t+d,t+d> dF_da;
        dF_da.template block<d,t+d>(0,0)=dd_da;
        dF_da.template block<t,t+d>(d,0)=-Cross_Product_Matrix(d_n)*dr_da+Cross_Product_Matrix(current_offset)*dd_da;
        return dF_da;
    }

    static Matrix<T,d,t+d> DD_DV(const RIGID_STRUCTURE<TV>& structure,const TV& object_offset,const TV& x1,const TV& x2,const TV& direction){
        T_SPIN a=structure.twist.angular;
        ROTATION<TV> orientation=structure.frame.orientation;
        orientation=ROTATION<TV>::From_Rotation_Vector(a).inverse()*orientation;
        TV offset=orientation*object_offset;
        auto norm_a=a.norm();
        auto dw_da=-sinc(norm_a/2)/4*a.transpose();
        auto a_norma=(norm_a>1e-8)?(TV)(a/norm_a):TV::UnitX();
        Matrix<T,3,3> dq_da=cos(norm_a/2)/2*a_norma*a_norma.transpose()+sinc(norm_a/2)/2*(Matrix<T,t,t>::Identity()-a_norma*a_norma.transpose());
        //std::cout<<dq_da<<std::endl;
        auto w=cos(norm_a/2);
        auto q=sinc(norm_a/2)*a/2;
        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);
        //dx_da.template block<d,t>(0,d)=Cross_Product_Matrix(offset).transpose();

        std::cout<<"DX_DA: "<<dx_da<<std::endl;
        
        auto distance=std::max((T)1e-3,(x2-x1).norm());
        auto dd_da=dx_da/distance-direction/cube(distance)*direction.transpose()*dx_da;
        std::cout<<"DD_DA: "<<dd_da<<std::endl;
        return dd_da;
    }

    static Matrix<T,1,t+d> DC_DA(const RIGID_STRUCTURE<TV>& structure,const TV& object_offset,const TV& x1,const TV& x2,const TV& direction){
        T_SPIN a=structure.twist.angular;
        ROTATION<TV> orientation=structure.frame.orientation;
        orientation=ROTATION<TV>::From_Rotation_Vector(a).inverse()*orientation;
        TV offset=orientation*object_offset;

        auto norm_a=a.norm();
        auto dw_da=-sinc(norm_a/2)/4*a.transpose();
        auto a_norma=(norm_a>1e-8)?(TV)(a/norm_a):TV::UnitX();
        Matrix<T,3,3> dq_da=cos(norm_a/2)/2*a_norma*a_norma.transpose()+sinc(norm_a/2)/2*(Matrix<T,t,t>::Identity()-a_norma*a_norma.transpose());
        //std::cout<<dq_da<<std::endl;
        auto w=cos(norm_a/2);
        auto q=sinc(norm_a/2)*a/2;
        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);
        //dx_da.template block<d,t>(0,d)=Cross_Product_Matrix(offset).transpose();

        //std::cout<<"DX_DA: "<<dx_da<<std::endl;
        
        auto distance=std::max((T)1e-3,(x2-x1).norm());
        //std::cout<<"Distance: "<<distance<<std::endl;
        //std::cout<<"Direction: "<<direction<<std::endl;
        TV normalized_direction=direction.normalized();
        //auto dd_da=dx_da/distance-direction/cube(distance)*direction.transpose()*dx_da;
        auto dd_da=dx_da/distance-normalized_direction/distance*normalized_direction.transpose()*dx_da;
        //std::cout<<"DD_DA: "<<dd_da<<std::endl;
        //std::cout<<normalized_direction<<std::endl;
        auto final=direction.transpose()*dd_da+normalized_direction.transpose()*(dx_da);
        //auto final=normalized_direction.transpose()*(dx_da);
        //std::cout<<"DC_DA: "<<final<<std::endl;
        return final;
    }

    static T Hessian(const T_SPIN& a1,const T_SPIN& a2,const TV& offset1,const TV& offset2,const TV& x1,const TV& x2,const TV& direction,int i1,int i2){
        static const T eps=1e-8;
        T scalar_zero(0);
        auto DX_DA=[](const T_SPIN& a,const TV& offset,int index){
            if(index<d){TV dx_da;dx_da.setZero();dx_da[index]=1;return dx_da;}
            T n_a=std::max(eps,a.norm());
            T c=cos(n_a/2);
            T s=sin(n_a/2);
            T dna_da=a[index-d]/n_a;
            T dw_da=-sinc(n_a/2)/4*dna_da;
            TV a_norma=(n_a>eps)?(TV)(a/n_a):TV::UnitX();
            TV dq_da=c/2*a_norma*dna_da+sinc(n_a/2)/2*(TV::Unit(index-d)-a_norma*dna_da);
            T w=c;
            TV q=sinc(n_a/2)*a/2;
            TV dx_da=-2*Cross_Product_Matrix(offset)*(dw_da*q+w*dq_da);
            return dx_da;
        };

        TV dr_da2=(i2<t+d?DX_DA(a1,offset1,i2):DX_DA(a2,offset2,i2-t-d));
        TV dr_da1=(i1<t+d?DX_DA(a1,offset1,i1):DX_DA(a2,offset2,i1-t-d));

        T distance=std::max(eps,(x2-x1).norm());
        TV dd_da1=dr_da1/distance-direction/cube(distance)*direction.transpose()*dr_da1;
        TV dd_da2=dr_da2/distance-direction/cube(distance)*direction.transpose()*dr_da2;

        TV r=x2-x1;
        T n_r=std::max(eps,r.norm());
        auto r_t=r.transpose();

        auto D2X_DAA=[](const T_SPIN& a,const TV& offset,int index){
            T n_a=std::max(eps,a.norm());
            T_SPIN da_da;da_da.setZero();da_da[index]=1;
            //T dna_da=(n_a>eps)?a[index]/n_a:scalar_zero;
            T dna_da=a[index]/n_a;
            //T da_na_da=(n_a>eps)?(1/n_a-a[index]/sqr(n_a)):scalar_zero;
            T da_na_da=1/n_a-a[index]/sqr(n_a);
            T c=cos(n_a/2);
            T s=sin(n_a/2);
            T dw_da=-sinc(n_a/2)/4*dna_da;
            TV a_norma=(n_a>eps)?(TV)(a/n_a):TV::UnitX();
            TV dq_da=c/2*a_norma*dna_da+sinc(n_a/2)/2*(da_da-a_norma*dna_da);
            T w=c;
            TV q=sinc(n_a/2)*a/2;

            T d2na_daa=1/n_a-1/sqr(n_a)*dna_da;
            T d2w_daa=-c/4*sqr(dna_da)-s/2*d2na_daa;
            TV d2q_daa=1/(2*n_a)*c*(2*da_da*dna_da-2*a_norma*sqr(dna_da)+a*d2na_daa)+1/n_a*s*(2*a/sqr(n_a)*sqr(dna_da)-2*da_da/n_a*dna_da-a_norma*d2na_daa-a/4*sqr(dna_da));
            auto oxstar=Cross_Product_Matrix(offset);
            TV d2x_daa=-2*oxstar*(d2w_daa*q+2*dw_da*dq_da+w*d2q_daa);
            return d2x_daa;
        };
        
        TV zero;zero.setZero();
        TV d2r_da1a2=((i1==i2&&i1>=(t+d+d))?D2X_DAA(a2,offset2,i1-t-d-d):zero)-((i1==i2&&i1<(t+d)&&i1>=d)?D2X_DAA(a1,offset1,i1-d):zero);
        TV d2d_da1a2=d2r_da1a2-dr_da1/cube(n_r)*r_t*dr_da2-r/cube(n_r)*(dr_da1.transpose()*dr_da2+r_t*d2r_da1a2)-dr_da2/cube(n_r)*r_t*dr_da1+3*r/std::pow(n_r,5)*r_t*dr_da2*r_t*dr_da1;
        TV d_n=direction.normalized();
        /*std::cout<<"d2r_da1a2: "<<d2r_da1a2.transpose()<<std::endl;
        std::cout<<"p1: "<<d2d_da1a2.transpose()*r<<std::endl;
        std::cout<<"p2: "<<dd_da1.transpose()*dr_da2<<std::endl;
        std::cout<<"p3: "<<dd_da2.transpose()*dr_da1<<std::endl;
        std::cout<<"p4: "<<d_n.transpose()*d2r_da1a2<<std::endl;*/
        T d2f_da1a2=(d2d_da1a2.transpose()*r+dd_da1.transpose()*dr_da2+dd_da2.transpose()*dr_da1+d_n.transpose()*d2r_da1a2).value();
        //std::cout<<"("<<i1<<","<<i2<<"): "<<d2f_da1a2<<std::endl;
        return d2f_da1a2;
    }

    static Matrix<T,2*(t+d),2*(t+d)> Full_Hessian(const T_SPIN& a1,const T_SPIN& a2,const TV& offset1,const TV& offset2,const TV& x1,const TV& x2,const TV& direction){
            // this is actually 2(t+d) x 2(t+d)
        Matrix<T,2*(t+d),2*(t+d)> full;
        for(int i=0;i<2*(t+d);i++){
            for(int j=0;j<2*(t+d);j++){
                //std::cout<<"Entry "<<i<<", "<<j<<std::endl;
                full(i,j)=Hessian(a1,a2,offset1,offset2,x1,x2,direction,i,j);
            }
        }
        return full;
    }
    
};
}
#endif
