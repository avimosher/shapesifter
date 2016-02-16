#ifndef __Relative_Position_Force__
#define __Relative_Position_Force__

#include <Math/RCF_NF.h>

namespace Mechanics{
template<class TV>
struct Relative_Position_Force
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;
    typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;

    template<int E,int ETYPE,int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void Apply(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms)
    {
        // do the F_NF version
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,Contract(F_NF<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),(VSIGN<E>::SIGN*error[E][ETYPE]).eval()),hessian_terms);
    }

    template<int E,int ETYPE,int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void Apply(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms)
    {
        // do the RCF_NF version
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,Contract(RCF_NF<TV,E>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),(VSIGN<E>::SIGN*error[E][ETYPE]).eval()),hessian_terms);
    }

    template<int E,int ETYPE,int V1,int VTYPE1>
    static void Apply_Part_One(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Apply,E,ETYPE,V1,VTYPE1),f,spin,offset,index,error,hessian_terms);
    }

    template<int E,int ETYPE>
    static void Derivatives(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Apply_Part_One,E,ETYPE),f,spin,offset,index,error,hessian_terms);
    }

    static void Build(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Derivatives),f,spin,offset,index,error,hessian_terms);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void Force_Velocity_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        TV term=(F_NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset))*error[E][ETYPE]*(T)VSIGN<E>::SIGN;
        Flatten_Matrix_Term<T,t+d,1,d,1>(index[V],constraint_index,VTYPE,0,term,force_constraint_terms);
        Flatten_Matrix_Term<T,1,t+d,1,d>(constraint_index,index[V],0,VTYPE,term.transpose(),constraint_force_terms);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void Force_Velocity_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        TV term=(RCF_NF<TV,E>::template First_Derivative<V,VTYPE>(f,spin,offset))*error[E][ETYPE]*(T)VSIGN<E>::SIGN;
        Flatten_Matrix_Term<T,t+d,1,d,1>(index[V],constraint_index,VTYPE,0,term,force_constraint_terms);
        Flatten_Matrix_Term<T,1,t+d,1,d>(constraint_index,index[V],0,VTYPE,term.transpose(),constraint_force_terms);
    }

    template<int E,int ETYPE>
    static void Force_Velocity_Derivatives(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Force_Velocity_Derivative,E,ETYPE),f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
    }

    static void Force_Velocity(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Force_Velocity_Derivatives),f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void First_Derivative(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[E],index[V],ETYPE,VTYPE,F_NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset).transpose()*scale*(T)VSIGN<E>::SIGN,force_terms);
    }

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void First_Derivative(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[E],index[V],ETYPE,VTYPE,RCF_NF<TV,E>::template First_Derivative<V,VTYPE>(f,spin,offset).transpose()*scale*(T)VSIGN<E>::SIGN,force_terms);
    }

    template<int E,int ETYPE>
    static void First_Derivatives_Second_Part(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(First_Derivative,E,ETYPE),index,scale,f,spin,offset,force_terms);
    }

    static void First_Derivatives(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(First_Derivatives_Second_Part),index,scale,f,spin,offset,force_terms);
    }

    template<int V1,int VTYPE1,int V2,int VTYPE2>
    static void Constraint_Second_Derivative(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,error*NF<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),hessian_terms);
    }

    template<int V1,int VTYPE1>
    static void Constraint_Second_Derivative_Part_Two(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Constraint_Second_Derivative,V1,VTYPE1),index,error,f,spin,offset,hessian_terms);
    }
    
    static void Constraint_Second_Derivatives(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Constraint_Second_Derivative_Part_Two),index,error,f,spin,offset,hessian_terms);
    }

    template<int V1,int VTYPE1,std::enable_if_t<VTYPE1==ANGULAR>* = nullptr>
    static void Evaluate(const std::array<int,2>& indices,const TV& f,const T scalar,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& right_hand_side){
        right_hand_side.template block<d,1>(indices[V1]*(d+d)+d,0)+=RCF_NF<TV,V1>::Evaluate(f,spin,offset)*(T)VSIGN<V1>::SIGN*scalar;
    }

    template<int V1,int VTYPE1,std::enable_if_t<VTYPE1==LINEAR>* = nullptr>
    static void Evaluate(const std::array<int,2>& indices,const TV& f,const T scalar,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& right_hand_side){
        right_hand_side.template block<d,1>(indices[V1]*(d+d),0)+=scalar*f/std::max((T)epsilon(),f.norm())*(T)VSIGN<V1>::SIGN;
    }

    static void Right_Hand_Sides(const std::array<int,2>& indices,const TV& f,const T scalar,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& right_hand_side){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Evaluate),indices,f,scalar,spin,offset,right_hand_side);
    }
};
}
#endif
