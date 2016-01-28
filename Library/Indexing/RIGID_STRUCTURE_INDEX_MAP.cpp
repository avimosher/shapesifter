#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_INDEX_MAP<TV>::
Compute_Constraint_Second_Derivatives(const Matrix<T,Dynamic,1>& force_balance_error,const std::array<int,2>& indices,int constraint_index,const T constraint_error,const T scalar_force,const TV& relative_position,const SparseMatrix<T>& f_scaling,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_constraint_terms,std::vector<Triplet<T>>& constraint_force_terms)
{
    std::array<T_SPIN,2> spins;
    std::array<TV,2> spun_offsets;
    for(int f=0,f_sign=-1;f<2;f++,f_sign+=2){
        for(int s1=0,s1_sign=-1;s1<2;s1++,s1_sign+=2){
            for(int s2=0,s2_sign=-1;s2<2;s2++,s2_sign+=2){
                // linear part WRT linear
                const T_TENSOR d2f_dv2=d2f_nf_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(relative_position,{s1_sign,s2_sign},spins,spun_offsets)*(T)f_sign*scalar_force;
                M_VxV f_times_d2f_dv2;f_times_d2f_dv2.setZero();
                for(int i=0;i<3;i++){
                    for(int j=0;j<3;j++){
                        for(int k=0;k<3;k++){
                            int f_index=indices[f]*(t+d)+k;
                            f_times_d2f_dv2(i,j)+=d2f_dv2(i,j,k)*f_scaling.coeff(f_index,f_index)*force_balance_error(f_index);}}}
                Flatten_Matrix_Term<T,t+d,t+d,d,d>(indices[s1],indices[s2],0,0,f_times_d2f_dv2,hessian_terms);}
            
            // second derivatives that have one force term and one velocity
            M_VxV df_dv=f_sign*df_nf_dVelocity<LINEARITY::LINEAR>(relative_position,s1_sign,spins[s1],spun_offsets[s1]);
            TV f_times_df_dv;f_times_df_dv.setZero();
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    int f_index=indices[f]*(t+d)+j;
                    f_times_df_dv(i)+=df_dv(j,i)*f_scaling.coeff(f_index,f_index)*force_balance_error(f_index);}}
            Flatten_Matrix_Term<T,t+d,1,d,1>(indices[s1],constraint_index,0,0,f_times_df_dv,force_constraint_terms);
            Flatten_Matrix_Term<T,1,t+d,1,d>(constraint_index,indices[s1],0,0,f_times_df_dv.transpose(),constraint_force_terms);}}
    // include second derivatives of the constraint equation itself (only double velocity)
    for(int s1=0,s1_sign=-1;s1<2;s1++,s1_sign+=2){
        for(int s2=0,s2_sign=-1;s2<2;s2++,s2_sign+=2){
            M_VxV full=d2n_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(relative_position,{s1_sign,s2_sign},{spins[s1],spins[s2]},{spun_offsets[s1],spun_offsets[s2]});
            Flatten_Matrix_Term<T,t+d,t+d,d,d>(indices[s1],indices[s2],0,0,constraint_error*full,hessian_terms);}}
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_INDEX_MAP)
