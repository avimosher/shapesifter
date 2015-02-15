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

        std::cout<<dx_da<<std::endl;
        
        auto distance=std::max((T)1e-3,(x2-x1).norm());
        auto dd_da=dx_da/distance-direction/cube(distance)*direction.transpose()*dx_da;
        std::cout<<dd_da<<std::endl;
        TV normalized_direction=direction.normalized();
        std::cout<<normalized_direction<<std::endl;
        auto final=normalized_direction.transpose()*dd_da+normalized_direction.transpose()*(dx_da);
        auto factor=normalized_direction.transpose()/distance-normalized_direction.transpose()*direction/cube(distance)*2*direction.transpose()+normalized_direction.transpose();
        std::cout<<"Factor: "<<factor<<std::endl;
        std::cout<<final<<std::endl;
        return final;
    }

    void Linearize(DATA<TV>& data,const T dt,const T time,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic);
    void Constraint_Satisfaction(DATA<TV>& data,const T dt,const T target_time,Matrix<T,Dynamic,1>& satisfaction);
    void Special(DATA<TV>& data,const T dt,const T time,SparseMatrix<T>& gradient);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("RELATIVE_POSITION_CONSTRAINT")
};
}
#endif
