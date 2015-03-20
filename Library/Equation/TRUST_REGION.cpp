#include <Data/DATA.h>
#include <Equation/EQUATION.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Equation/TRUST_REGION.h>
#include <Force/FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/MATH.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> TRUST_REGION<TV>::
TRUST_REGION()
{
    function_scale_factor=1;
    f*=function_scale_factor;
    gk*=function_scale_factor;
    norm_gk=gk.norm();
    prec=1e-8;
    contract_factor=.5;
    expand_factor=3;
    contract_threshold=.2;
    expand_threshold_ap=.8;
    expand_threshold_rad=.8;
    trust_iterations=200;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Resize_Vectors()
{
    nvars=equation->System_Size();
    xk.setZero(nvars);
    gk.resize(nvars);
    sk.resize(nvars);
    yk.resize(nvars);
    try_x.resize(nvars);
    try_g.resize(nvars);

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
Linearize(SIMULATION<TV>& simulation,const T dt,const T time)
{
    DATA<TV>& data=simulation.data;
    FORCE<TV>& force=simulation.force;

    // zero velocities
    current_velocities.resize(data.Velocity_DOF(),1);current_velocities.setZero();
    data.Unpack_Velocities(current_velocities);

    // zero forces
    force.Pack_Forces(solve_forces);
    solve_forces.setZero();
    force.Unpack_Forces(solve_forces);

    // store positions
    data.Pack_Positions(positions);

    equation->Linearize(data,force,dt,time,true);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Linearize_Around()
{
    DATA<TV>& data=stored_simulation->data;
    FORCE<TV>& force=stored_simulation->force;
    // sk is solve_vector
    Vector solve_velocities=sk.block(0,0,data.Velocity_DOF(),1);
    solve_forces.Set(sk.block(data.Velocity_DOF(),0,sk.rows()-data.Velocity_DOF(),1));
    data.Unpack_Positions(positions);
    force.Increment_Forces(solve_forces,1);
    data.Unpack_Velocities(current_velocities+solve_velocities);
    data.Step();
    equation->Linearize(data,force,dt,time,false);
    force.Increment_Forces(solve_forces,-1);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Increment_X()
{
    DATA<TV>& data=stored_simulation->data;
    FORCE<TV>& force=stored_simulation->force;
    // sk is solve_vector
    Vector solve_velocities=sk.block(0,0,data.Velocity_DOF(),1);
    solve_forces.Set(sk.block(data.Velocity_DOF(),0,sk.rows()-data.Velocity_DOF(),1));
    force.Increment_Forces(solve_forces,1);
    current_velocities+=solve_velocities;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    // based on trustOptim implementation
    int iteration=0;
    auto status=CONTINUE;

    stored_simulation=&simulation;
    this->dt=dt;
    this->time=time;

    iteration=0;
    max_iterations=100;
    radius=1;
    min_radius=1e-9;
    preconditioner_refresh_frequency=1;

    Linearize(simulation,dt,time);
    Get_F(xk,f);
    Get_Gradient(xk,gk);
    Resize_Vectors();
    norm_gk=gk.norm();
    Update_Hessian();
    PrecondLLt.analyzePattern(Bk);
    Update_Preconditioner();
    do{
        iteration++;
        status=Update_One_Step();
        std::cout<<"norm_gk: "<<norm_gk<<" norm_gk/sqrt(nvars): "<<norm_gk/sqrt(T(nvars))<<std::endl;
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
            simulation.Write("Frame "+std::to_string(simulation.current_frame)+" substep "+std::to_string(iteration));
        }
        if(status==CONTRACT){status=CONTINUE;}
    }while(status==CONTINUE);
    std::cout<<"SOLVE STEPS: "<<iteration<<std::endl;
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
    
    /*for(int j=0;j<BB.outerSize();j++){
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
    else{alpha=beta/2;}*/
    int ii=0;
    do{
        ii++;
        PrecondLLt.factorize(BB);
        if(PrecondLLt.info()==Eigen::Success){
            success=true;
        }
        /*else{
            alpha=std::max(2*alpha,beta/2)-alpha;
            for(int j=0;j<nvars;j++){
                BB.coeffRef(j,j)+=alpha;
            }
            }*/
    }while(!success);

    //PrecondLLt.analyzePattern(BB);
    //PrecondLLt.factorize(BB);

    /*Matrix<T,Dynamic,Dynamic> L(PrecondLLt.matrixL());
    Matrix<T,Dynamic,Dynamic> Lt(PrecondLLt.matrixU());
    std::cout<<"LLt"<<std::endl;
    std::cout<<PrecondLLt.permutationPinv()*L*Lt*PrecondLLt.permutationP()<<std::endl;
    std::cout<<"BB"<<std::endl;
    std::cout<<BB<<std::endl;*/
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
        //try_x=xk+sk;
        Linearize_Around();
        Get_F(try_x,try_f);
        std::cout<<"Old f: "<<f<<" try f: "<<try_f<<std::endl;
        if(finite(try_f)){
            try_f*=function_scale_factor;
            ared=f-try_f;
            gs=gk.dot(sk);
            sBs=sk.dot(Bk.template selfadjointView<Lower>()*sk);
            pred=-(gs+sBs/2);
            if(pred<0){step_status=ENEGMOVE;}
            ap=ared/pred;
            std::cout<<"AP: "<<ap<<" ared: "<<ared<<" pred: "<<pred<<" radius: "<<radius<<std::endl;
        }
        else{step_status=FAILEDCG;}
    }
    if(step_status!=FAILEDCG && step_status!=ENEGMOVE){
        if(ap>contract_threshold){
            Get_Gradient(try_x,try_g);
            if(finite(try_g.norm())){
                try_g*=function_scale_factor;
                yk=try_g-gk;
                f=try_f;
                Increment_X();
                //xk+=sk;
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
            step_status=CONTRACT;
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
    std::cout<<"CG reason: "<<CG_stop_reason<<" iterations: "<<num_CG_iterations<<std::endl;
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
    //equation->Linearize_Around(x);
    f=equation->Evaluate();
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Get_Gradient(const Vector& x,Vector& g)
{
    //equation->Linearize_Around(x);
    g=equation->Gradient();
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
UPz(const Preconditioner& X,const Vector& v,Vector& out)
{
    //out=X.permutationP()*v;
    //out=X.matrixU().template triangularView<Upper>()*out;
    // it looks like X.matrixU() is L...
    out=X.permutationP()*v;
    out=X.matrixU()*out;
    out=X.matrixL()*out;
    out=X.permutationPinv()*out;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar TRUST_REGION<TV>::
Find_Tau(const Vector& z,const Vector& d)
{
    if(d[0]>.7974&&d[0]<.7975){
        SparseMatrix<T> L=PrecondLLt.matrixL();
        SparseMatrix<T> Lt=PrecondLLt.matrixU();
        std::cout<<"LLt"<<std::endl;
        std::cout<<L*Lt<<std::endl;
        std::cout<<"Lt"<<std::endl;
        std::cout<<Lt.template triangularView<Upper>()<<std::endl;
        std::cout<<"m_matrix"<<std::endl;
        //std::cout<<PrecondLLt.rawMatrix()<<std::endl;
        SparseMatrix<T> BB(nvars,nvars);
        BB=Bk.template selfadjointView<Lower>();
        std::cout<<"BB"<<std::endl;
        std::cout<<BB<<std::endl;
    }
    UPz(PrecondLLt,d,wd);
    UPz(PrecondLLt,z,wz);
    
    /*
    T d2=wd.squaredNorm();
    T z2=z.squaredNorm();
    T zd=wd.dot(wz);
    
    T root=zd*zd-d2*(z2-radius*radius);
    T tau=(sqrt(root)-zd)/d2;
    */

    T pCd=z.dot(wd);
    T dCd=d.dot(wd);
    T pCp=z.dot(wz);
    T tau=(-2*pCd+sqrt(pCd*pCd-4*dCd*(pCp-radius*radius)))/(2*dCd);
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
DEFINE_AND_REGISTER_PARSER(TRUST_REGION,void)
{
    auto step=std::make_shared<TRUST_REGION<TV>>();
    step->equation=new NONLINEAR_EQUATION<TV>();
    simulation.evolution.push_back(step);
    return 0;
}
