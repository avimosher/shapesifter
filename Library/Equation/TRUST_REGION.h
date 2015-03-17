#ifndef __TRUST_REGION__
#define __TRUST_REGION__

#include <Evolution/EVOLUTION_STEP.h>

namespace Mechanics{
    template<class TV> class EQUATION;
    
template<class TV>
class TRUST_REGION:public EVOLUTION_STEP<TV>
{
    typedef typename TV::Scalar T;
    typedef Matrix<T,Dynamic,1> Vector;
    typedef SimplicialLLT<SparseMatrix<T>> Preconditioner;
public:
    enum STATUS{UNKNOWN,CONTINUE,SUCCESS,EMAXITER,ETOLG,MOVED,EXPAND,CONTRACT,FAILEDCG,ENEGMOVE};
    EQUATION<TV>* equation;
    SparseMatrix<T> Bk;
    Preconditioner PrecondLLt;
    Vector xk;
    Vector gk;
    Vector yk;
    Vector sk;
    Vector try_x;
    Vector try_g;
    Vector zj;
    Vector rj;
    Vector dj;
    Vector zj_old;
    Vector yj;
    Vector wd;
    Vector wz;
    T try_f,gs;
    T f;
    T norm_gk,ared,pred,sBs,ap,norm_sk_scaled;
    
    T radius;
    T min_radius;
    T tol;
    T prec;
    T contract_factor;
    T expand_factor;
    T contract_threshold;
    T expand_threshold_ap;
    T expand_threshold_rad;
    T function_scale_factor;
    int preconditioner_refresh_frequency;
    int max_iterations,num_CG_iterations,trust_iterations;
    int nvars;
    std::string CG_stop_reason;

    TRUST_REGION();
    ~TRUST_REGION(){}

    void Step(SIMULATION<TV>& simulation,const T dt,const T time);
    void Resize_Vectors();
    void Linearize(SIMULATION<TV>& simulation,const T dt,const T time);
    void Update_Preconditioner();
    void Update_Hessian();
    STATUS Update_One_Step();
    void Solve_Trust_CG(Vector& pk);
    T Get_Norm_Sk(const Preconditioner& X);
    void Get_F(const Vector& x,T& f);
    void Get_Gradient(const Vector& x,Vector& g);
    void UPz(const Preconditioner& X,const Vector& v,Vector& out);
    T Find_Tau(const Vector& z,const Vector& d);
    void Get_FDF(const Vector& x,T& f,Vector& g);
    DEFINE_TYPE_NAME("TRUST_REGION")
};
}
#endif
