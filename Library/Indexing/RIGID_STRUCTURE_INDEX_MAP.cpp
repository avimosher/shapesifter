#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Math/Relative_Position_Force.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_INDEX_MAP<TV>::
Compute_Constraint_Second_Derivatives(const Matrix<T,Dynamic,1>& force_balance_error,const std::array<int,2>& indices,int constraint_index,const T constraint_error,const T scalar_force,const TV& relative_position,const std::array<T_SPIN,2>& spin,const std::array<TV,2>& offset,const SparseMatrix<T>& f_scaling,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms)
{
    std::array<std::array<TV,2>,2> error;
    std::array<std::array<TV,2>,2> unforced_error;
    for(int s1=0;s1<2;s1++){
        for(int i=0;i<3;i++){
            int f_index=indices[s1]*(t+d)+i;
            unforced_error[s1][0][i]=force_balance_error(f_index)*f_scaling.coeff(f_index,f_index);
            error[s1][0][i]=scalar_force*unforced_error[s1][0][i];
            int r_index=indices[s1]*(t+d)+i+3;
            unforced_error[s1][1][i]=force_balance_error(r_index)*f_scaling.coeff(r_index,r_index);
            error[s1][1][i]=scalar_force*unforced_error[s1][1][i];}}
    Relative_Position_Force<TV>::Build(relative_position,spin,offset,indices,error,hessian_terms);
    Relative_Position_Force<TV>::Force_Velocity(relative_position,spin,offset,indices,constraint_index,unforced_error,force_constraint_terms,constraint_force_terms);
    Relative_Position_Force<TV>::Constraint_Second_Derivatives(indices,constraint_error,relative_position,spin,offset,hessian_terms);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_INDEX_MAP)
