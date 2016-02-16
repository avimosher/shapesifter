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

    //*********************************************
    // d2F/dv1/dv2

    template<int E,int ETYPE,int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void Apply(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms)
    {
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,Contract(F_NF<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),(VSIGN<E>::SIGN*error[E][ETYPE]).eval()),hessian_terms);}

    template<int E,int ETYPE,int V1,int VTYPE1,int V2,int VTYPE2,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void Apply(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms)
    {
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,Contract(RCF_NF<TV,E>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),(VSIGN<E>::SIGN*error[E][ETYPE]).eval()),hessian_terms);}

    template<int E,int ETYPE,int V1,int VTYPE1>
    static void Apply_Part_One(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Apply,E,ETYPE,V1,VTYPE1),f,spin,offset,index,error,hessian_terms);}

    template<int E,int ETYPE>
    static void Derivatives(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Apply_Part_One,E,ETYPE),f,spin,offset,index,error,hessian_terms);}

    static void Second_Derivatives(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Derivatives),f,spin,offset,index,error,hessian_terms);}

    //*********************************************
    // d2F/dv/df

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void Force_Velocity_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        TV term=(F_NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset))*error[E][ETYPE]*(T)VSIGN<E>::SIGN;
        Flatten_Matrix_Term<T,t+d,1,d,1>(index[V],constraint_index,VTYPE,0,term,force_constraint_terms);
        Flatten_Matrix_Term<T,1,t+d,1,d>(constraint_index,index[V],0,VTYPE,term.transpose(),constraint_force_terms);}

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void Force_Velocity_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        TV term=(RCF_NF<TV,E>::template First_Derivative<V,VTYPE>(f,spin,offset))*error[E][ETYPE]*(T)VSIGN<E>::SIGN;
        Flatten_Matrix_Term<T,t+d,1,d,1>(index[V],constraint_index,VTYPE,0,term,force_constraint_terms);
        Flatten_Matrix_Term<T,1,t+d,1,d>(constraint_index,index[V],0,VTYPE,term.transpose(),constraint_force_terms);}

    template<int E,int ETYPE>
    static void Force_Velocity_Derivatives(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Force_Velocity_Derivative,E,ETYPE),f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);}

    static void Force_Velocity(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const std::array<int,2>& index,const int constraint_index,const std::array<std::array<TV,2>,2>& error,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Force_Velocity_Derivatives),f,spin,offset,index,constraint_index,error,force_constraint_terms,constraint_force_terms);}

    //*********************************************
    // dF/dv

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static void Balance_Velocity_Derivative(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[E],index[V],ETYPE,VTYPE,F_NF<TV>::template First_Derivative<V,VTYPE>(f,spin,offset).transpose()*scale*(T)VSIGN<E>::SIGN,force_terms);}

    template<int E,int ETYPE,int V,int VTYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static void Balance_Velocity_Derivative(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[E],index[V],ETYPE,VTYPE,RCF_NF<TV,E>::template First_Derivative<V,VTYPE>(f,spin,offset).transpose()*scale*(T)VSIGN<E>::SIGN,force_terms);}

    template<int E,int ETYPE>
    static void Balance_Velocity_Derivatives_Second_Part(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Balance_Velocity_Derivative,E,ETYPE),index,scale,f,spin,offset,force_terms);}

    static void Balance_Velocity_Derivatives(const std::array<int,2>& index,const T scale,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& force_terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Balance_Velocity_Derivatives_Second_Part),index,scale,f,spin,offset,force_terms);}

    //*********************************************
    // dF/df

    template<int E,int ETYPE,std::enable_if_t<ETYPE==LINEAR>* = nullptr>
    static TV Balance_Force_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return f.normalized();
    }

    template<int E,int ETYPE,std::enable_if_t<ETYPE==ANGULAR>* = nullptr>
    static TV Balance_Force_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return (ROTATION<TV>::From_Rotation_Vector(spin[E])*offset[E]).cross(f.normalized());
    }

    template<int E,int ETYPE>
    static void Store_Balance_Force_Derivative(const std::array<int,2>& index,int constraint_index,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& terms){
        Flatten_Matrix_Term<T,t+d,1,d,1>(index[E],constraint_index,ETYPE,0,Balance_Force_Derivative<E,ETYPE>(f,spin,offset)*(T)VSIGN<E>::SIGN,terms);
    }

    // for a particular constraint.  There should be four derivatives.
    static void Balance_Force_Derivatives(const std::array<int,2>& index,int constraint_index,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Store_Balance_Force_Derivative),index,constraint_index,f,spin,offset,terms);
    }


    //*********************************************
    // F

    template<int V1,int VTYPE1,std::enable_if_t<VTYPE1==ANGULAR>* = nullptr>
    static void Evaluate(const std::array<int,2>& indices,const TV& f,const T scalar,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& right_hand_side){
        right_hand_side.template block<d,1>(indices[V1]*(d+d)+d,0)+=RCF_NF<TV,V1>::Evaluate(f,spin,offset)*(T)VSIGN<V1>::SIGN*scalar;}

    template<int V1,int VTYPE1,std::enable_if_t<VTYPE1==LINEAR>* = nullptr>
    static void Evaluate(const std::array<int,2>& indices,const TV& f,const T scalar,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& right_hand_side){
        right_hand_side.template block<d,1>(indices[V1]*(d+d),0)+=scalar*f/std::max((T)epsilon(),f.norm())*(T)VSIGN<V1>::SIGN;}

    static void Right_Hand_Sides(const std::array<int,2>& indices,const TV& f,const T scalar,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,Matrix<T,Dynamic,1>& right_hand_side){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Evaluate),indices,f,scalar,spin,offset,right_hand_side);}
};
}
#endif
