#ifndef __Relative_Position_Constraint_Term__
#define __Relative_Position_Constraint_Term__

#include <Math/NF.h>

namespace Mechanics{
template<class TV>
struct Relative_Position_Constraint_Term
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum {d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,d,d> M_VxV;

    static T Evaluate(const TV& relative_position,const T target_distance){
        return relative_position.norm()-target_distance;
    }

    static void Constraint_Velocity_Derivatives(){
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
