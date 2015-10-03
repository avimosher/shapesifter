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
    //typedef SimplicialLLT<SparseMatrix<T>> Preconditioner;
    typedef IncompleteCholesky<T> Preconditioner;
public:
    enum STATUS{UNKNOWN,CONTINUE,SUCCESS,EMAXITER,ETOLG,MOVED,EXPAND,CONTRACT,FAILEDCG,ENEGMOVE,NEGRATIO};
    EQUATION<TV>* equation;
    SparseMatrix<T> hessian;
    Preconditioner preconditioner;
    Vector xk,gk,sk,try_g,zj,rj,dj,zj_old,yj,wd,wz;
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
    Vector positions;
    SIMULATION<TV>* stored_simulation;
    T dt;
    T time;

    TRUST_REGION();
    ~TRUST_REGION(){}

    void Step(SIMULATION<TV>& simulation,const T dt,const T time);
    void Resize_Vectors();
    void Linearize(SIMULATION<TV>& simulation,const T dt,const T time);
    void Linearize_Around();
    void Increment_X();;
    void Update_Preconditioner();
    void Update_Hessian();
    STATUS Update_One_Step();
    void Solve_Trust_CG(Vector& pk);
    void Multiply(const Preconditioner& X,const Vector& v,Vector& out);
    T Find_Tau(const Vector& z,const Vector& d);
    DEFINE_TYPE_NAME("TRUST_REGION")
};
}
#endif
