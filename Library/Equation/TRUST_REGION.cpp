#include <Equation/EQUATION.h>
#include <Equation/TRUST_REGION.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/MATH.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> TRUST_REGION<TV>::
TRUST_REGION()
{
    Get_FDF(xk,f,gk);
    f*=function_scale_factor;
    gk*=function_scale_factor;
    norm_gk=gk.norm();
    zj.setZero(nvars);
    rj.setZero(nvars);
    dj.setZero(nvars);
    zj_old.setZero(nvars);
    yj.setZero(nvars);
    wd.resize(nvars);
    wz.resize(nvars);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    // based on trustOptim implementation
    int iteration=0;
    auto status=CONTINUE;
    
    do{
        iteration++;
        status=Update_One_Step();
        if(norm_gk/sqrt(T(nvars))<=prec){
            status=SUCCESS;
        }
        if(iteration>=max_iterations){
            status=EMAXITER;
        }
        if(radius<=min_radius){ // trust region collapse
            status=ETOLG;
        }

        // update Hessian
        if(status==MOVED || status==EXPAND){
            Update_Hessian();
            if(iteration%preconditioner_refresh_frequency==0){
                Update_Preconditioner();
            }
            status=CONTINUE;
        }
        if(status==CONTRACT){status=CONTINUE;}
    }while(status==CONTINUE);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Update_Preconditioner()
{
    bool success=false;
    T alpha,beta,bmin;

    Vector TT(nvars);
    SparseMatrix<T> BB(nvars,nvars);
    BB=Bk.template selfadjointView<Lower>();
    
    for(int j=0;j<BB.outerSize();j++){
        TT(j)=sqrt(BB.innerVector(j).dot(BB.innerVector(j)));
    }

    for(int j=0;j<BB.outerSize();j++){
        for(typename SparseMatrix<T>::InnerIterator it(BB,j);it;++it){
            BB.coeffRef(it.row(),j)*=1/sqrt(TT(it.row())*TT(j));
        }
    }

    beta=sqrt(BB.cwiseAbs2().sum());
    bmin=BB.coeff(0,0);
    for(int j=0;j<nvars;j++){
        bmin=std::min(bmin,BB.coeff(j,j));
    }
    
    if(bmin>0){alpha=0;}
    else{alpha=beta/2;}

    int ii=0;
    do{
        ii++;
        PrecondLLt.factorize(BB);
        if(PrecondLLt.info()==Eigen::Success){
            success=true;
        }
        else{
            alpha=std::max(2*alpha,beta/2)-alpha;
            for(int j=0;j<nvars;j++){
                BB.coeffRef(j,j)+=alpha;
            }
        }
    }while(!success);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Update_Hessian()
{
    // call to NONLINEAR_EQUATION
    Bk=equation->Hessian();
    Bk*=function_scale_factor;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TRUST_REGION<TV>::STATUS TRUST_REGION<TV>::
Update_One_Step()
{
    auto step_status=UNKNOWN;
    Solve_Trust_CG(sk);
    norm_sk_scaled=Get_Norm_Sk(PrecondLLt);
    if(!finite(norm_sk_scaled)){step_status=FAILEDCG;}
    else{
        try_x=xk+sk;
        Get_F(try_x,try_f);
        if(finite(try_f)){
            try_f*=function_scale_factor;
            ared=f-try_f;
            gs=gk.dot(sk);
            sBs=sk.dot(Bk.template selfadjointView<Lower>()*sk);
            pred=-(gs+sBs/2);
            if(pred<0){step_status=ENEGMOVE;}
            ap=ared/pred;
        }
        else{step_status=FAILEDCG;}
    }
    if(step_status!=FAILEDCG && step_status!=ENEGMOVE){
        if(ap>contract_threshold){
            Gradient(try_x,try_g);
            if(finite(try_g.norm())){
                try_g*=function_scale_factor;
                yk=try_g-gk;
                f=try_f;
                xk+=sk;
                gk=try_g;
                norm_gk=gk.norm();
                if(ap>expand_threshold_ap && norm_sk_scaled>=expand_threshold_rad*radius){
                    step_status=EXPAND;
                }
                else{step_status=MOVED;}
            }
            else{step_status=FAILEDCG;}
        }
        else{
            step_status=CONTRACT;
        }
    }

    switch(step_status){
        case CONTRACT:
        case FAILEDCG:
        case ENEGMOVE:
            radius*=contract_factor;
            break;
        case EXPAND:
            radius*=expand_factor;
            break;
    };
    return step_status;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Solve_Trust_CG(Vector& pk)
{
    T norm_rj,dot_ry,dot_ry_old,norm_zj,aj,bj,tau,dBd,norm_gk;
    int j;
    T crit;
    
    zj.setZero();
    rj=-gk;
    UPz(PrecondLLt,rj,wd);
    norm_rj=wd.norm();
    UPz(PrecondLLt,gk,wd);
    norm_gk=wd.norm();

    // Solve LL'y=r
    yj=PrecondLLt.solve(rj);
    dj=yj;
    
    std::stringstream reason;
    for(j=0;j<trust_iterations;j++){
        dBd=dj.dot(Bk.template selfadjointView<Lower>()*dj);
        if(dBd<=0){
            tau=Find_Tau(zj,dj);
            pk=zj+tau*dj;
            num_CG_iterations=j+1;
            reason<<"Negative curvature";
            break;}

        aj=rj.dot(yj)/dBd;
        zj_old=zj;
        zj.noalias()+=aj*dj;

        UPz(PrecondLLt,zj,wd);
        norm_zj=wd.norm();

        if(norm_zj>=radius){
            // find tau>=0 s.t. p intersects trust region
            tau=Find_Tau(zj_old,dj);
            pk=zj_old+tau*dj;
            num_CG_iterations=j+1;
            reason<<"Intersect TR bound";
            break;}

        dot_ry=rj.dot(yj);
        rj.noalias()-=aj*(Bk.template selfadjointView<Lower>()*dj).eval();
        UPz(PrecondLLt,rj,wd);
        norm_rj=wd.norm();
        crit=norm_rj/norm_gk;
        
        if(crit<tol){
            pk=zj;
            num_CG_iterations=j+1;
            reason<<"Reached tolerance";
            break;}

        dot_ry_old=dot_ry;
        
        //updating yj
        yj=PrecondLLt.solve(rj);
        dot_ry=rj.dot(yj);
        bj=dot_ry/dot_ry_old;
        dj*=bj;
        dj.noalias()+=yj;
    }
    
    if(j>=trust_iterations){
        pk=zj;
        num_CG_iterations=j;
        reason<<"Exceeded max CG iterations";
    }

    CG_stop_reason=reason.str();
    return;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar TRUST_REGION<TV>::
Get_Norm_Sk(const Preconditioner& X)
{
    T res=(X.matrixU()*(X.permutationPinv()*sk).eval()).norm();
    return res;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Get_F(const Vector& x,T& f)
{
    equation->Linearize_Around(x);
    f=equation->Evaluate();
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Gradient(const Vector& x,Vector& g)
{
    equation->Linearize_Around(x);
    g=equation->Gradient();
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
UPz(const Preconditioner& X,const Vector& v,Vector& out)
{
    out=X.permutationP()*v;
    out=X.matrixU().template triangularView<Upper>()*out;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar TRUST_REGION<TV>::
Find_Tau(const Vector& z,const Vector& d)
{
    UPz(PrecondLLt,d,wd);
    UPz(PrecondLLt,z,wz);
    
    T d2=wd.squaredNorm();
    T z2=z.squaredNorm();
    T zd=wd.dot(wz);
    
    T root=zd*zd-d2*(z2-radius*radius);
    T tau=(sqrt(root)-zd)/d2;
    return tau;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Get_FDF(const Vector& x,T& f,Vector& g)
{
    // TODO: CALL MODEL
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(TRUST_REGION)
