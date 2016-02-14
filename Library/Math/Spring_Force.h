#ifndef __Spring_Force__
#define __Spring_Force__

#include <Math/RCF_NF.h>

namespace Mechanics{
template<class TV>
struct Spring_Force
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;

    template<int E,int ETYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static TV Evaluate(const T k,const T target,const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        TV f=F<TV>::Evaluate(x,spin,offset);
        T nf=std::max((T)1e-8,f.norm());
        return -k*(nf-target)*f/nf*VSIGN<E>::SIGN;
    }

    template<int E,int ETYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static T_SPIN Evaluate(const T k,const T target,const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        TV f=F<TV>::Evaluate(x,spin,offset);
        T nf=std::max((T)1e-8,f.norm());
        return (ROTATION<TV>::From_Rotation_Vector(spin[E])*offset[E]).cross(-k*(nf-target)*f/nf*VSIGN<E>::SIGN);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static M_VxV First_Derivative(const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=std::max((T)1e-8,f.norm());
        return -k*VSIGN<E>::SIGN*(NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset)*F_NF<TV>::Evaluate(f,spin,offset).transpose()+(nf-target)*F_NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset));
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static M_VxV First_Derivative(const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        T nf=std::max((T)1e-8,f.norm());
        return -k*VSIGN<E>::SIGN*(RCF_NF<TV,E>::template First_Derivative<V,VTYPE>(f,spin,offset)*(nf-target)+NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset)*RCF_NF<TV,E>::Evaluate(f,spin,offset).transpose());
    }

    template<int E,int ETYPE,int V,int VTYPE>
    static void Matrix_Term(const std::array<int,2>& index,const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[E],index[V],ETYPE,VTYPE,First_Derivative<E,ETYPE,V,VTYPE>(k,target,f,spin,offset).transpose(),force_terms);
    }

    template<int E,int ETYPE>
    static void First_Derivatives(const std::array<int,2>& index,const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Matrix_Term<E,ETYPE,0,LINEAR>(index,k,target,f,spin,offset,force_terms);
        Matrix_Term<E,ETYPE,0,ANGULAR>(index,k,target,f,spin,offset,force_terms);
        Matrix_Term<E,ETYPE,1,LINEAR>(index,k,target,f,spin,offset,force_terms);
        Matrix_Term<E,ETYPE,1,ANGULAR>(index,k,target,f,spin,offset,force_terms);
    }

    static void Derivatives(const std::array<int,2>& index,const T k,const T target,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        First_Derivatives<0,LINEAR>(index,k,target,f,spin,offset,force_terms);
        First_Derivatives<0,ANGULAR>(index,k,target,f,spin,offset,force_terms);
        First_Derivatives<1,LINEAR>(index,k,target,f,spin,offset,force_terms);
        First_Derivatives<1,ANGULAR>(index,k,target,f,spin,offset,force_terms);
    }

    template<int E,int ETYPE>
    static void Evaluate_Term(const std::array<int,2>& index,const T k,const T target,const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& error){
        error.template block<d,1>(index[E]*(t+d)+ETYPE*d,0)+=Evaluate<E,ETYPE>(k,target,x,spin,offset);
    }

    static void Evaluate(const std::array<int,2>& index,const T k,const T target,const std::array<TV,2>& x,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& error){
        Evaluate_Term<0,LINEAR>(index,k,target,x,spin,offset,error);
        Evaluate_Term<0,ANGULAR>(index,k,target,x,spin,offset,error);
        Evaluate_Term<1,LINEAR>(index,k,target,x,spin,offset,error);
        Evaluate_Term<1,ANGULAR>(index,k,target,x,spin,offset,error);
    }
};
}
#endif
