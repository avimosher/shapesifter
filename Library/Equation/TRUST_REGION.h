#ifndef __TRUST_REGION__
#define __TRUST_REGION__

#include <Evolution/EVOLUTION_STEP.h>
#include <Force/STORED_FORCE.h>
#include <Eigen/IterativeSolvers>

namespace Mechanics{
template<class TV> class EQUATION;
    
template<class TV>
class TRUST_REGION:public EVOLUTION_STEP<TV>
{
    typedef typename TV::Scalar T;
    typedef Matrix<T,Dynamic,1> Vector;
    typedef IncompleteCholesky<T> Preconditioner;
public:
    enum STATUS{UNKNOWN,CONTINUE,SUCCESS,EMAXITER,ETOLG,MOVED,EXPAND,CONTRACT,FAILEDCG,ENEGMOVE,NEGRATIO};
    EQUATION<TV>* equation;
    SparseMatrix<T> hessian;
    SparseMatrix<T> jacobian;
    Preconditioner preconditioner;
    Vector gk,sk,try_g,zj,rj,dj,zj_old,yj,wd,wz,inverse_scale,Jsk,componentwise_prediction,rhs,try_rhs;
    T f;
    T norm_gk;
    
    T radius;
    T min_radius;
    T tol;
    T precision;
    T contract_factor;
    T expand_factor;
    T contract_threshold;
    T expand_threshold_ap;
    T expand_threshold_rad;
    int preconditioner_refresh_frequency;
    int max_iterations,num_CG_iterations,trust_iterations;
    int nvars;
    std::string CG_stop_reason;

    STORED_FORCE<T> solve_forces;
    Vector current_velocities;
    Vector candidate_velocities;
    Vector positions;

    TRUST_REGION();
    ~TRUST_REGION(){}

    void Step(SIMULATION<TV>& simulation,const T dt,const T time);
    void Linearize(SIMULATION<TV>& simulation,const T dt,const T time);
    void Linearize_Around(SIMULATION<TV>& simulation,const T dt,const T time);
    void Increment_X(SIMULATION<TV>& simulation);
    void Update_Preconditioner(bool identity);
    void Update_Hessian();
    STATUS Update_One_Step(SIMULATION<TV>& simulation,const T dt,const T time);
    void Solve_Trust_Conjugate_Gradient(Vector& pk);
    void Multiply(const Preconditioner& X,const Vector& v,Vector& out);
    T Norm(const Preconditioner& X,const Vector& v,Vector& scratch);
    T Find_Tau(const Vector& z,const Vector& d);
    void Check_Derivative(SIMULATION<TV>& simulation,const T dt,const T time);
    DEFINE_TYPE_NAME("TRUST_REGION")
};
}
#endif
