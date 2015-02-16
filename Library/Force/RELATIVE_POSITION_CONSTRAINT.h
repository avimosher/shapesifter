#ifndef __RELATIVE_POSITION_CONSTRAINT__
#define __RELATIVE_POSITION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

namespace Mechanics{

template<class TV>
class RELATIVE_POSITION_CONSTRAINT : public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
public:
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    using FORCE_TYPE<TV>::stored_forces;
    struct CONSTRAINT{
        int s1;
        TV v1;
        int s2;
        TV v2;
        T target_distance;

        template<class Archive>
        void serialize(Archive& archive)
        {archive(s1,v1,s2,v2,target_distance);}
    };
    std::vector<CONSTRAINT> constraints;

    RELATIVE_POSITION_CONSTRAINT(){}
    ~RELATIVE_POSITION_CONSTRAINT(){}

    template<class Archive>
    void serialize(Archive& archive)
    {archive(constraints);}

    Matrix<T,1,t+d> DC_DA(const T_SPIN& a,const TV& offset,const TV& x1,const TV& x2,const TV& direction){
        auto norm_a=a.norm();
        auto dw_da=-sinc(norm_a/2)/4*a.transpose();
        auto a_norma=(norm_a>1e-8)?(TV)(a/norm_a):TV::UnitX();
        Matrix<T,3,3> dq_da=cos(norm_a/2)/2*a_norma*a_norma.transpose()+sinc(norm_a/2)/2*(Matrix<T,t,t>::Identity()-a_norma*a_norma.transpose());
        std::cout<<dq_da<<std::endl;
        auto w=cos(norm_a/2);
        auto q=sinc(norm_a/2)*a/2;
        Matrix<T,d,t+d> dx_da; // identity for translation parts, zero for rotation
        dx_da.template block<d,d>(0,0).setIdentity();
        dx_da.template block<d,t>(0,d)=-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);
        //dx_da.template block<d,t>(0,d)=Cross_Product_Matrix(offset).transpose();

        std::cout<<dx_da<<std::endl;
        
        auto distance=std::max((T)1e-3,(x2-x1).norm());
        auto dd_da=dx_da/distance-direction/cube(distance)*direction.transpose()*dx_da;
        std::cout<<dd_da<<std::endl;
        TV normalized_direction=direction.normalized();
        std::cout<<normalized_direction<<std::endl;
        auto final=direction.transpose()*dd_da+normalized_direction.transpose()*(dx_da);
        //auto final=normalized_direction.transpose()*(dx_da);
        std::cout<<final<<std::endl;
        return final;
    }


    void Hessian(const T_SPIN& a,const TV& offset,const TV& x1,const TV& x2,const TV& direction,int index){
        auto n_a=a.norm();
        T_SPIN da_da;da_da.setZero();da_da[index]=1;
        auto dna_da=a[index]/n_a; // TODO: index
        auto c=cos(n_a/2);
        auto s=sin(n_a/2);
        auto dw_da=-sinc(n_a/2)/4*dna_da;
        auto a_norma=(n_a>1e-8)?(TV)(a/n_a):TV::UnitX();
        Matrix<T,3,1> dq_da=c/2*a_norma*dna_da+sinc(n_a/2)/2*(da_da-a_norma/n_a*dna_da);
        auto w=c;
        auto q=sinc(n_a/2)*a/2;
        /*Matrix<T,d,t+d> dr_da2; // TODO: dr_da actually needs to select which one we're using, sign, etc
        dr_da2.template block<d,d>(0,0).setIdentity();
        dr_da2.template block<d,t>(0,d)=-2*Cross_Product_Matrix(offset)*(q*dw_da+w*dq_da);*/
        Matrix<T,d,1> dr_da2;

        auto dr_da1=dr_da2;

        auto distance=std::max((T)1e-3,(x2-x1).norm());
        auto dd_da1=dr_da1/distance-direction/cube(distance)*direction.transpose()*dr_da1;
        auto dd_da2=dr_da2/distance-direction/cube(distance)*direction.transpose()*dr_da2;
        

        auto r=x2-x1;
        auto n_r=r.norm();
        auto r_t=r.transpose();
        auto d2na_daa=1/n_a-1/sqr(n_a)*dna_da;
        auto d2w_daa=-c/4*sqr(dna_da)-s/2*d2na_daa;
        auto d2q_daa=1/(2*n_a)*c*(2*da_da*dna_da-2*a/n_a*sqr(dna_da)+a*d2na_daa)+1/n_a*s*(2*a/sqr(n_a)*sqr(dna_da)-2*da_da/n_a*dna_da-a/n_a*d2na_daa-a/4*sqr(dna_da));
        auto oxstar=Cross_Product_Matrix(offset);
        auto d2x_daa=-2*oxstar*(d2w_daa*q+2*dw_da*dq_da+w*d2q_daa);
        auto d2r_da1a2=d2x_daa; // TODO: deal with using the right x, a
        auto p1=d2r_da1a2;
        auto p2=dr_da1/cube(n_r)*r_t*dr_da2;
        auto p3=r/cube(n_r)*(dr_da1.transpose()*dr_da2+r_t*d2r_da1a2);
        auto p4=dr_da2/cube(n_r)*r_t*dr_da1;
        auto p5=3*r/std::pow(n_r,5)*r_t*dr_da2*r_t*dr_da1;
        auto d2d_da1a2=d2r_da1a2-dr_da1/cube(n_r)*r_t*dr_da2-r/cube(n_r)*(dr_da1.transpose()*dr_da2+r_t*d2r_da1a2)-dr_da2/cube(n_r)*r_t*dr_da1+3*r/std::pow(n_r,5)*r_t*dr_da2*r_t*dr_da1;
        auto d=direction.normalized();
        auto d2f_da1a2=d2d_da1a2.transpose()*r+dd_da1.transpose()*dr_da2+dd_da2.transpose()*dr_da1+d.transpose()*d2r_da1a2;


    }

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    void Constraint_Satisfaction(DATA<TV>& data,const T dt,const T target_time,Matrix<T,Dynamic,1>& satisfaction);
    void Special(DATA<TV>& data,const T dt,const T time,SparseMatrix<T>& gradient);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("RELATIVE_POSITION_CONSTRAINT")
};
}
#endif
