#include <Data/DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/EQUATION.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Equation/TRUST_REGION.h>
#include <Force/FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/LOG.h>
#include <Utilities/MATH.h>
#include <Utilities/RANDOM.h>
#include <Eigen/Eigenvalues>
#include <iomanip>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> TRUST_REGION<TV>::
TRUST_REGION()
{
    precision=1e-6;
    contract_factor=.25;
    expand_factor=2.5;
    contract_threshold=.4;
    expand_threshold=.8;
    expand_threshold_rad=.8;
    trust_iterations=20;
    tol=1e-8;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Linearize(SIMULATION<TV>& simulation,const T dt,const T time)
{
    DATA<TV>& data=simulation.data;
    FORCE<TV>& force=simulation.force;

    // zero velocities
    equation->Initialize(data,force);
    current_velocities.resize(equation->Velocity_DOF(),1);current_velocities.setZero();
    equation->Unpack_Velocities(data,current_velocities);

    // zero forces
    force.Pack_Forces(solve_forces);
    solve_forces.setZero();
    force.Unpack_Forces(solve_forces);

    // store positions
    data.Pack_Positions(positions);

    equation->Linearize(data,force,dt,time,true);
    equation->RHS(rhs);
    force.Pack_Forces(solve_forces);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Linearize_Around(SIMULATION<TV>& simulation,const T dt,const T time,const Vector& solve_vector)
{
    DATA<TV>& data=simulation.data;
    FORCE<TV>& force=simulation.force;
    // sk is solve_vector
    int velocity_dof=equation->Velocity_DOF();
    Vector solve_velocities=solve_vector.block(0,0,velocity_dof,1);
    solve_forces.Set(solve_vector.block(velocity_dof,0,solve_vector.rows()-velocity_dof,1));
    data.Unpack_Positions(positions);
    force.Increment_Forces(solve_forces,1);
    
    candidate_velocities=current_velocities+solve_velocities;
    equation->Unpack_Velocities(data,candidate_velocities);
    data.Step();
    equation->Linearize(data,force,dt,time,false);
    force.Increment_Forces(solve_forces,-1);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Increment_X(SIMULATION<TV>& simulation)
{
    // sk is solve_vector
    int velocity_dof=equation->Velocity_DOF();
    Vector solve_velocities=sk.block(0,0,velocity_dof,1);
    solve_forces.Set(sk.block(velocity_dof,0,sk.rows()-velocity_dof,1));
    simulation.force.Increment_Forces(solve_forces,1);
    simulation.force.Pack_Forces(solve_forces);
    current_velocities+=solve_velocities;
    rhs=try_rhs;
    
    // store errors
    T one_over_maxabs=1;
    if(rhs.rows()>0){
        T maxabs=rhs.array().abs().maxCoeff();
        if(!maxabs){maxabs=1;}
        one_over_maxabs=1/maxabs;}
    equation->Store_Errors(simulation.data,rhs.block(0,0,velocity_dof,1)*one_over_maxabs);
    solve_forces.Set(rhs.block(velocity_dof,0,rhs.rows()-velocity_dof,1)*one_over_maxabs);
    simulation.force.Store_Errors(solve_forces);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    static int call_count=0;
    static int total_steps=0;
    // loosely based on trustOptim implementation
    status=CONTINUE;

    iteration=0;
    max_iterations=100;
    radius=1;
    min_radius=1e-9;
    preconditioner_refresh_frequency=1;

    Linearize(simulation,dt,time);
    f=equation->Evaluate();
    equation->Gradient(gk);
    norm_gk=gk.norm();
    Update_Hessian(false);
    Update_Preconditioner(false);

    static int failed_radius=0;
    do{
        iteration++;
        LOG::cout<<std::endl<<"BEGINNING STEP "<<iteration<<std::endl;
        status=Update_One_Step(simulation,dt,time);
        LOG::cout<<"norm_gk: "<<norm_gk<<" norm_gk/sqrt(nvars): "<<norm_gk/sqrt(T(nvars))<<std::endl;
        //if(norm_gk/sqrt(T(nvars))<=precision && f<=precision){status=SUCCESS;}
        if(norm_gk/sqrt(T(nvars))<=precision){status=SUCCESS;}
        if(iteration>=max_iterations){status=EMAXITER;}
        if(radius<=min_radius){ // trust region collapse
            status=ETOLG;
            failed_radius++;
        }

        // update Hessian
        if(status==MOVED || status==EXPAND){
            //Update_Hessian(status!=EXPAND);
            Update_Hessian(false);
            //Check_Derivative(simulation,dt,time);
            if(simulation.force.Equations_Changed() || iteration%preconditioner_refresh_frequency==0){Update_Preconditioner(false);}
            status=CONTINUE;}
        if(simulation.substeps){
            std::string frame_name="Frame "+std::to_string(simulation.current_frame)+" substep "+std::to_string(iteration)+" real "+std::to_string(int(status==CONTINUE))+ " f "+std::to_string(f);
            simulation.Write(frame_name);}
        if(status==CONTRACT){status=CONTINUE;}
    }while(status==CONTINUE);
    std::string frame_name="End frame "+std::to_string(simulation.current_frame)+" substep "+std::to_string(iteration)+" real "+std::to_string(int(status==CONTINUE))+ " f "+std::to_string(f);
    simulation.Write(frame_name);
    //Check_Derivative(simulation,dt,time);
    LOG::cout<<"SOLVE STEPS: "<<iteration<<" Failed due to radius: "<<failed_radius<<std::endl;
    call_count++;
    total_steps+=iteration;
    LOG::cout<<"Current average: "<<total_steps/(T)call_count<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Update_Preconditioner(bool identity)
{
    /*Matrix<T,Dynamic,Dynamic> dense(hessian);
    EigenSolver<Matrix<T,Dynamic,Dynamic>> es(dense,false);
    LOG::cout<<es.eigenvalues()<<std::endl;*/
    if(!identity){
        preconditioner.compute(hessian);}
    if(identity || preconditioner.info()!=ComputationInfo::Success){
        LOG::cout<<"Preconditioner computation failed; using identity"<<std::endl;
        SparseMatrix<T> BB(nvars,nvars);
        BB.setIdentity();
        preconditioner.compute(BB);}
    /*inverse_scale.resize(nvars);
    for(int i=0;i<nvars;i++){
    inverse_scale(i)=1/preconditioner.scalingS()(i);}*/
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Update_Hessian(bool use_accurate_hessian)
{
    if(use_accurate_hessian){
        equation->Accurate_Hessian(hessian);}
    else{
        equation->Hessian(hessian);}
    equation->Jacobian(jacobian);
    nvars=hessian.rows();
}
///////////////////////////////////////////////////////////////////////
template<class M>
void Print_No_Angular(const M& m)
{
    for(int i=0;i<m.rows();i++){
        if(i%6>2){continue;}
        for(int j=0;j<m.cols();j++){
            if(j%6>2){continue;}
            LOG::cout<<std::setw(13)<<m(i,j)<<" ";
        }
        LOG::cout<<std::endl;
    }
}
template<class TV> typename TRUST_REGION<TV>::STATUS TRUST_REGION<TV>::
Update_One_Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    auto step_status=UNKNOWN;
    T try_f,step_quality,predicted_reduction;
    /*LOG::cout<<"Jacobian adjoint: "<<std::endl;
    Print_No_Angular(Matrix<T,Dynamic,Dynamic>(jacobian.adjoint()));
    LOG::cout<<"Error: "<<std::endl;
    LOG::cout<<"Un-adjusted hessian sk:"<<std::endl;*/
    //hessian=jacobian.adjoint()*jacobian;
    //Solve_Trust_Conjugate_Gradient(sk);
    //equation->Hessian(hessian);
    Solve_Trust_Conjugate_Gradient(sk);
    //Solve_Trust_MINRES(sk);
    T norm_sk_scaled=Norm(preconditioner,sk,wd);
    if(!std::isfinite(norm_sk_scaled)){step_status=FAILEDCG;}
    else{
        Linearize_Around(simulation,dt,time,sk);
        try_f=equation->Evaluate();
        if(std::isfinite(try_f)){
            T actual_reduction=f-try_f;
            T gs=gk.dot(sk);
            T sBs=sk.dot(hessian.template selfadjointView<Lower>()*sk);
            predicted_reduction=-(gs+sBs/2);
            //Jsk=jacobian*sk;
            //componentwise_prediction=gk.cwiseProduct(Jsk)+Jsk.cwiseProduct(Jsk)/2;
            equation->RHS(try_rhs);
            //LOG::cout<<"Try RHS: "<<std::endl<<try_rhs<<std::endl;
            //LOG::cout<<"Componentwise quality: "<<componentwise_prediction.sum()<<std::endl<<rhs.Diff(try_rhs).cwiseQuotient(componentwise_prediction)<<std::endl;
            //int index;
            //LOG::cout<<"Min value at index "<<index<<" is "<<componentwise_prediction.minCoeff(&index)<<std::endl;
            //LOG::cout<<"Max value at index "<<index<<" is "<<componentwise_prediction.maxCoeff(&index)<<std::endl;
            if(predicted_reduction<0){step_status=ENEGMOVE;}
            step_quality=actual_reduction/predicted_reduction;
            LOG::cout<<"Candidate error with value "<<try_f<<":"<<std::endl;
            /*Matrix<T,Dynamic,1> error;
            equation->RHS(error);
            LOG::cout<<"Current error: "<<std::endl;
            Print_No_Angular(error.transpose());*/
            LOG::cout<<"AP: "<<step_quality<<" old f: "<<f<<" try f: "<<try_f<<" ared: "<<actual_reduction<<" pred: "<<predicted_reduction<<" radius: "<<radius<<" gs: "<<gs<<" sBs: "<<sBs<<" norm_sk_scaled: "<<norm_sk_scaled<<std::endl;}
        else{step_status=FAILEDCG;}}
    if(step_status!=FAILEDCG && step_status!=ENEGMOVE){
        if(step_quality>contract_threshold){
            equation->Gradient(try_g);
            if(std::isfinite(try_g.norm())){
                f=try_f;
                Increment_X(simulation);
                gk=try_g;
                norm_gk=gk.norm();
                /*LOG::cout<<"Resolved error:"<<std::endl;
                Matrix<T,Dynamic,1> error;
                equation->RHS(error);
                Print_No_Angular(error.transpose());*/
                if(step_quality>expand_threshold){// && norm_sk_scaled>=expand_threshold_rad*radius){
                    step_status=EXPAND;}
                else{step_status=MOVED;}}
            else{step_status=FAILEDCG;}}
        else if(step_quality<0){step_status=NEGRATIO;}
        else{step_status=CONTRACT;}}

    LOG::cout<<"Step status: "<<step_status<<std::endl;
    switch(step_status){
        case NEGRATIO:{
            T gksk=gk.dot(sk);
            T gamma_bad=(1-contract_threshold)*gksk/((1-contract_threshold)*(f+gksk+contract_threshold*(f-predicted_reduction)-try_f));
            radius=std::min(contract_factor*norm_sk_scaled,std::max((T)0.0625,gamma_bad)*radius);
            step_status=CONTRACT;
            break;}
        case CONTRACT:
        case FAILEDCG:
        case ENEGMOVE:
            step_status=CONTRACT;
            //radius=norm_sk_scaled*contract_factor;
            radius*=contract_factor;
            break;
        case EXPAND:
            //radius=std::max(expand_factor*norm_sk_scaled,radius);
            radius*=expand_factor;
            break;
        default:
            break;
    };
    return step_status;
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar TRUST_REGION<TV>::
Norm(const Preconditioner& preconditioner,const Vector& v,Vector& scratch)
{
    if(preconditioner.permutationP().rows() == v.rows()){
        scratch=preconditioner.permutationP()*v;}
    else{scratch=v;}
    scratch=preconditioner.matrixL().adjoint().template triangularView<Upper>()*scratch;
    return scratch.norm();
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Solve_Trust_MINRES(Vector& sol)
{
    std::array<T,2> alpha,beta,res;
    std::array<Vector,2> v,d;
    int n=hessian.rows();
    for(int i=0;i<2;i++){
        v[i].resize(n);
        v[i].setZero();
        d[i].resize(n);
        d[i].setZero();
    }
    T norm_A=0;
    T cond_A=1;
    T c=-1,s=0; // givens

    T gamma_min=1e99;
    std::array<T,2> delta1,delta2,ep,gamma1,gamma2;
    delta1[1]=0;
    Vector pk,tk,old_sol;
    sol.resize(n);sol.setZero();
    T tau=0;
    T epsilon=1e-8;
    beta[0]=0;
    T norm_rhs=gk.norm();
    beta[1]=norm_rhs;
    v[1]=-gk/norm_rhs;
    res[0]=norm_rhs;
    tau=norm_rhs;

    auto sign=[&](T x) {return (std::fabs(x)<epsilon?0:x/std::fabs(x));};

    int i;
    std::stringstream reason;
    for(i=0;i<trust_iterations;i++){
        int cur=(i+1)%2,next=i%2;

        pk=hessian*v[cur];
        alpha[cur]=v[cur].dot(pk);
        pk-=alpha[cur]*v[cur];
        v[next]=pk-beta[cur]*v[next];
        beta[next]=v[next].norm();
        if(fabs(beta[next])>epsilon){
            v[next]/=beta[next];}

        delta2[cur]=c*delta1[cur]+s*alpha[cur];
        gamma1[cur]=s*delta1[cur]-c*alpha[cur];

        ep[next]=s*beta[next];
        delta1[next]=-c*beta[next];


        T a=gamma1[cur],b=beta[next];
        if(fabs(b)<epsilon){
            s=0;
            gamma2[cur]=fabs(a);
            if(fabs(a)<epsilon){
                c=1;}
            else{
                c=sign(a);}
        }
        else if(fabs(a)<epsilon){
            c=0;
            s=sign(b);
            gamma2[cur]=fabs(b);
        }
        else if(fabs(b)>fabs(a)){
            T t=a/b;
            s=sign(b)/sqrt(1+sqr(t));
            c=s*t;
            gamma2[cur]=b/s;}
        else{
            T t=b/a;
            c=sign(a)/sqrt(1+sqr(t));
            s=c*t;
            gamma2[cur]=a/c;}

        tau=c*res[next];
        res[cur]=s*res[next];

        if(i==0){norm_A=sqrt(sqr(alpha[cur])+sqr(beta[next]));}
        else{
            T tnorm=sqrt(sqr(alpha[cur])+sqr(beta[next])+sqr(beta[cur]));
            norm_A=std::max(norm_A,tnorm);}

        if(fabs(gamma2[cur])>epsilon){
            d[cur]=(v[cur]-delta2[cur]*d[next]-ep[cur]*d[cur])/gamma2[cur];

            old_sol=sol;
            sol+=tau*d[cur];
            if(sol.norm()>=radius){
                LOG::cout<<"tau was originally "<<tau<<std::endl;
                LOG::cout<<"Direction "<<d[cur].transpose()<<std::endl;
                d[cur]*=tau;
                tau=Find_Tau(old_sol,d[cur]);
                LOG::cout<<"Chosen tau is "<<tau<<std::endl;
                sol=old_sol+tau*d[cur];
                reason<<"Intersect TR bound";
                break;
            }
            gamma_min=std::min(gamma_min,gamma2[cur]);
            cond_A=norm_A/gamma_min;
        }

        LOG::cout<<"residual: "<<res[i%2]<<std::endl;
        if(res[i%2]/norm_rhs<tol){
            reason<<"Reached tolerance";
            break;
        }
    }
    CG_stop_reason=reason.str();
    LOG::cout<<"MINRES reason: "<<CG_stop_reason<<" iterations: "<<i<<std::endl;
    LOG::cout<<"solutn: "<<sol.transpose()<<std::endl;
    LOG::cout<<"result: "<<(hessian*sol).transpose()<<std::endl;
    LOG::cout<<"actual: "<<-gk.transpose()<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Solve_Trust_Conjugate_Gradient(Vector& pk)
{
    T dot_ry,dot_ry_old,aj,tau,dBd,p_norm_gk;
    int j;
    zj.resize(hessian.rows());
    zj.setZero();
    rj=-gk;
    p_norm_gk=Norm(preconditioner,gk,wd);
    T local_tol=std::min((T).5,(T)sqrt(p_norm_gk))*p_norm_gk;
    LOG::cout<<"Local tol: "<<local_tol<<" tol: "<<tol<<std::endl;

    // Solve LL'y=r
    yj=preconditioner.solve(rj);
    dj=yj;
    
    std::stringstream reason;
    for(j=0;j<trust_iterations;j++){
        dBd=dj.dot(hessian.template selfadjointView<Lower>()*dj);
        if(dBd<=0){
            tau=Find_Tau(zj,dj);
            pk.noalias()=zj+tau*dj;
            num_CG_iterations=j+1;
            reason<<"Negative curvature: "<<dBd;
            break;}

        aj=rj.dot(yj)/dBd;
        zj_old=zj;
        zj.noalias()+=aj*dj;

        if(Norm(preconditioner,zj,wd)>=radius){
            // find tau>=0 s.t. p intersects trust region
            tau=Find_Tau(zj_old,dj);
            pk.noalias()=zj_old+tau*dj;
            num_CG_iterations=j+1;
            reason<<"Intersect TR bound";
            break;}

        dot_ry=rj.dot(yj);
        rj.noalias()-=aj*(hessian.template selfadjointView<Lower>()*dj).eval();
        
        if(Norm(preconditioner,rj,wd)/p_norm_gk<local_tol){
            pk=zj;
            num_CG_iterations=j+1;
            reason<<"Reached tolerance";
            break;}

        dot_ry_old=dot_ry;
        
        //updating yj
        yj=preconditioner.solve(rj);
        dot_ry=rj.dot(yj);
        dj*=dot_ry/dot_ry_old;
        dj.noalias()+=yj;}
    
    if(j>=trust_iterations){
        pk=zj;
        num_CG_iterations=j;
        reason<<"Exceeded max CG iterations";}

    CG_stop_reason=reason.str();
    LOG::cout<<"CG reason: "<<CG_stop_reason<<" iterations: "<<num_CG_iterations<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Multiply(const Preconditioner& X,const Vector& v,Vector& out)
{
    if(X.permutationP().rows() == v.rows()){
        out=X.permutationP()*v;}
    else{out=v;}
    out=X.matrixL().adjoint().template triangularView<Upper>()*out;
    out=X.matrixL().template triangularView<Lower>()*out;
    if(X.permutationP().rows() == v.rows()){
        out=X.permutationP().inverse()*out;}
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar TRUST_REGION<TV>::
Find_Tau(const Vector& z,const Vector& d)
{
    Multiply(preconditioner,d,wd);
    Multiply(preconditioner,z,wz);
    
    T pCd=z.dot(wd);
    T dCd=d.dot(wd);
    T pCp=z.dot(wz);
    return (-2*pCd+sqrt(4*pCd*pCd-4*dCd*(pCp-radius*radius)))/(2*dCd);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Check_Derivative(SIMULATION<TV>& simulation,const T dt,const T time)
{
    DATA<TV>& data=simulation.data;
    FORCE<TV>& force=simulation.force;

    T epsilon=1e-2;
    Matrix<T,Dynamic,1> variables=equation->Get_Unknowns(data,force);
    variables.setZero();
    if(!sk.rows()){sk.resize(variables.rows());sk.setZero();}
    Linearize_Around(simulation,dt,time,variables);

    data.random.Direction(variables);
    T f0=equation->Evaluate();
    Matrix<T,Dynamic,1> gradient;equation->Gradient(gradient);
    SparseMatrix<T> h;equation->Hessian(h);
    auto Evaluate_Step_Error = [&](T eps){
        Linearize_Around(simulation,dt,time,eps*variables);

        T f1=equation->Evaluate();
        T predicted_delta_f=gradient.dot(eps*variables)+(T).5*eps*eps*variables.transpose()*h*variables;
        T error=f1-f0-predicted_delta_f;
        LOG::cout<<"Error for "<<eps<<": "<<error<<std::endl;
        return error;
    };
    T last_error=Evaluate_Step_Error(epsilon);
    int divisors=8;
    for(int i=0;i<divisors;i++){
        epsilon/=2;
        T new_error=Evaluate_Step_Error(epsilon);
        LOG::cout<<"Ratio: "<<last_error/new_error<<std::endl;
        last_error=new_error;
    }
    //LOG::cout<<"Ratio: "<<Evaluate_Step_Error(epsilon)/Evaluate_Step_Error(epsilon/2)<<std::endl;
    Linearize_Around(simulation,dt,time,sk);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(TRUST_REGION)
DEFINE_AND_REGISTER_PARSER(TRUST_REGION,void)
{
    auto step=std::make_shared<TRUST_REGION<TV>>();
    step->equation=new NONLINEAR_EQUATION<TV>();
    Parse_String(node["name"],step->name);
    Parse_Scalar(node["precision"],step->precision,step->precision);
    Parse_Scalar(node["contract_threshold"],step->contract_threshold,step->contract_threshold);
    Parse_Scalar(node["expand_threshold"],step->expand_threshold,step->expand_threshold);
    Parse_Scalar(node["contract_factor"],step->contract_factor,step->contract_factor);
    Parse_Scalar(node["expand_factor"],step->expand_factor,step->expand_factor);
    Parse_Scalar(node["trust_iterations"],step->trust_iterations,step->trust_iterations);
    Parse_Scalar(node["tol"],step->tol,step->tol);
    simulation.evolution.push_back(step);
    return 0;
}
