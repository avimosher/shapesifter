#ifndef __Relative_Position_Constraint_Term__
#define __Relative_Position_Constraint_Term__

#include <Math/NF.h>
#include <Math/RXO.h>

namespace Mechanics{
template<class TV>
struct Relative_Position_Constraint_Term
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,1,d> TV_T;
    typedef Matrix<T,d,d> M_VxV;

    //*********************************************
    // C

    static T Evaluate(const TV& relative_position,const T target_distance){
        return relative_position.norm()-target_distance;
    }

    //*********************************************
    // dC/dv

    template<int V,int VTYPE,std::enable_if_t<VTYPE==LINEAR>* = nullptr>
    static TV_T Constraint_Velocity_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return clamped_normalize(f).transpose();
    }

    template<int V,int VTYPE,std::enable_if_t<VTYPE==ANGULAR>* = nullptr>
    static TV_T Constraint_Velocity_Derivative(const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset){
        return (RXO<TV,V>::template First_Derivative<V,VTYPE>(f,spin,offset)*clamped_normalize(f)).transpose();
    }

    template<int V,int VTYPE>
    static void Store_Constraint_Velocity_Derivative(const std::array<int,2>& index,int constraint_index,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& terms){
        Flatten_Matrix_Term<T,1,t+d,1,d>(constraint_index,index[V],0,VTYPE,Constraint_Velocity_Derivative<V,VTYPE>(f,spin,offset)*(T)VSIGN<V>::SIGN,terms);
    }

    // for a particular constraint.  There should be four derivatives.
    static void Constraint_Velocity_Derivatives(const std::array<int,2>& index,int constraint_index,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Store_Constraint_Velocity_Derivative),index,constraint_index,f,spin,offset,terms);
    }
    

    //*********************************************
    // d2C/dv1/dv2

    template<int V1,int VTYPE1,int V2,int VTYPE2>
    static void Constraint_Second_Derivative(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        Flatten_Matrix_Term<T,t+d,t+d,d,d>(index[V1],index[V2],VTYPE1,VTYPE2,error*NF<TV>::template Second_Derivative<V1,VTYPE1,V2,VTYPE2>(f,spin,offset),hessian_terms);}

    template<int V1,int VTYPE1>
    static void Constraint_Second_Derivative_Part_Two(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(EXPAND_FUNCTION(Constraint_Second_Derivative,V1,VTYPE1),index,error,f,spin,offset,hessian_terms);}
    
    static void Constraint_Second_Derivatives(const std::array<int,2>& index,const T error,const TV& f,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,std::vector<Triplet<T>>& hessian_terms){
        EXPAND_BODIES_TYPES(WRAP_FUNCTION(Constraint_Second_Derivative_Part_Two),index,error,f,spin,offset,hessian_terms);}

};
}
#endif
