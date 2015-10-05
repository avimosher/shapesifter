#include <Data/DATA.h>
#include <Equation/EQUATION.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Equation/TRUST_REGION.h>
#include <Force/FORCE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/LOG.h>
#include <Utilities/MATH.h>
#include <Eigen/Eigenvalues>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> TRUST_REGION<TV>::
TRUST_REGION()
{
    precision=1e-6;
    contract_factor=.25;
    expand_factor=2.5;
    contract_threshold=.4;
    expand_threshold_ap=.8;
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
    current_velocities.resize(data.Velocity_DOF(),1);current_velocities.setZero();
    data.Unpack_Velocities(current_velocities);

    // zero forces
    force.Pack_Forces(solve_forces);
    solve_forces.setZero();
    force.Unpack_Forces(solve_forces);

    // store positions
    data.Pack_Positions(positions);

    equation->Linearize(data,force,dt,time,true);
    force.Pack_Forces(solve_forces);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Linearize_Around(SIMULATION<TV>& simulation,const T dt,const T time)
{
    DATA<TV>& data=simulation.data;
    FORCE<TV>& force=simulation.force;
    // sk is solve_vector
    int velocity_dof=data.Velocity_DOF();
    Vector solve_velocities=sk.block(0,0,velocity_dof,1);
    solve_forces.Set(sk.block(velocity_dof,0,sk.rows()-velocity_dof,1));
    data.Unpack_Positions(positions);
    force.Increment_Forces(solve_forces,1);
    data.Unpack_Velocities(current_velocities+solve_velocities);
    data.Step();
    equation->Linearize(data,force,dt,time,false);
    force.Increment_Forces(solve_forces,-1);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Increment_X(SIMULATION<TV>& simulation)
{
    // sk is solve_vector
    int velocity_dof=simulation.data.Velocity_DOF();
    Vector solve_velocities=sk.block(0,0,velocity_dof,1);
    solve_forces.Set(sk.block(velocity_dof,0,sk.rows()-velocity_dof,1));
    simulation.force.Increment_Forces(solve_forces,1);
    simulation.force.Pack_Forces(solve_forces);
    current_velocities+=solve_velocities;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    // loosely based on trustOptim implementation
    int iteration=0;
    auto status=CONTINUE;

    iteration=0;
    max_iterations=100;
    radius=1;
    min_radius=1e-9;
    preconditioner_refresh_frequency=1;

    Linearize(simulation,dt,time);
    f=equation->Evaluate();
    equation->Gradient(gk);
    norm_gk=gk.norm();
    Update_Hessian();
    Update_Preconditioner();

    static int failed_radius=0;
    do{
        iteration++;
        status=Update_One_Step(simulation,dt,time);
        LOG::cout<<"norm_gk: "<<norm_gk<<" norm_gk/sqrt(nvars): "<<norm_gk/sqrt(T(nvars))<<std::endl;
        if(norm_gk/sqrt(T(nvars))<=precision && f<=precision){status=SUCCESS;}
        if(iteration>=max_iterations){status=EMAXITER;}
        if(radius<=min_radius){ // trust region collapse
            status=ETOLG;
            failed_radius++;}

        // update Hessian
        if(status==MOVED || status==EXPAND){
            Update_Hessian();
            if(iteration%preconditioner_refresh_frequency==0){Update_Preconditioner();}
            status=CONTINUE;}
        if(simulation.substeps){
            std::string frame_name="Frame "+std::to_string(simulation.current_frame)+" substep "+std::to_string(iteration)+" radius "+std::to_string(radius)+" real "+std::to_string(int(status==CONTINUE))+ " f "+std::to_string(f);
            std::cout<<"WRITING FRAME "<<frame_name<<std::endl;
            simulation.Write(frame_name);}
        if(status==CONTRACT){status=CONTINUE;}
    }while(status==CONTINUE);
    LOG::cout<<"SOLVE STEPS: "<<iteration<<" Failed due to radius: "<<failed_radius<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Update_Preconditioner()
{
    nvars=hessian.rows();
#if 0
    Vector TT(nvars);
    SparseMatrix<T> BB(nvars,nvars);
    BB.setIdentity();
    preconditioner.compute(BB);
#else
    preconditioner.compute(hessian);
    if(preconditioner.info()!=ComputationInfo::Success){
        LOG::cout<<"Preconditioner computation failed; using identity"<<std::endl;
        SparseMatrix<T> BB(nvars,nvars);
        BB.setIdentity();
        preconditioner.compute(BB);
    }
    /*
    Matrix<T,Dynamic,Dynamic> dense(hessian);
    EigenSolver<Matrix<T,Dynamic,Dynamic>> es(dense);
    LOG::cout<<"Matrix: "<<std::endl<<dense<<std::endl;
    LOG::cout<<"Eigenvalues: "<<es.eigenvalues()<<std::endl;
    Matrix<T,Dynamic,1> first,second;
    first=Matrix<T,Dynamic,1>::Constant(nvars,1,1);
    LOG::cout<<preconditioner.solve(first).transpose()<<std::endl;
    LOG::cout<<dense.jacobiSvd(ComputeThinU | ComputeThinV).solve(first).transpose()<<std::endl;*/
#endif
}
///////////////////////////////////////////////////////////////////////
template<class TV> void TRUST_REGION<TV>::
Update_Hessian()
{
    equation->Hessian(hessian);
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TRUST_REGION<TV>::STATUS TRUST_REGION<TV>::
Update_One_Step(SIMULATION<TV>& simulation,const T dt,const T time)
{
    auto step_status=UNKNOWN;
    T try_f,step_quality,predicted_reduction;
    Solve_Trust_Conjugate_Gradient(sk);
    T norm_sk_scaled=preconditioner._norm(sk,wd);
    //LOG::cout<<"sk: norm: "<<norm_sk_scaled<<std::endl<<sk.transpose()<<std::endl;
    if(!finite(norm_sk_scaled)){step_status=FAILEDCG;}
    else{
        Linearize_Around(simulation,dt,time);
        try_f=equation->Evaluate();
        if(finite(try_f)){
            T actual_reduction=f-try_f;
            T gs=gk.dot(sk);
            T sBs=sk.dot(hessian.template selfadjointView<Lower>()*sk);
            predicted_reduction=-(gs+sBs/2);
            if(predicted_reduction<0){step_status=ENEGMOVE;}
            step_quality=actual_reduction/predicted_reduction;
            LOG::cout<<"AP: "<<step_quality<<" old f: "<<f<<" try f: "<<try_f<<" ared: "<<actual_reduction<<" pred: "<<predicted_reduction<<" radius: "<<radius<<" gs: "<<gs<<" sBs: "<<sBs<<" norm_sk_scaled: "<<norm_sk_scaled<<std::endl;}
        else{step_status=FAILEDCG;}}
    if(step_status!=FAILEDCG && step_status!=ENEGMOVE){
        if(step_quality>contract_threshold){
            equation->Gradient(try_g);
            if(finite(try_g.norm())){
                f=try_f;
                Increment_X(simulation);
                gk=try_g;
                norm_gk=gk.norm();
                //tol=std::min(0.5,sqrt(norm_gk))*norm_gk;
                if(step_quality>expand_threshold_ap && norm_sk_scaled>=expand_threshold_rad*radius){
                    step_status=EXPAND;}
                else{step_status=MOVED;}}
            else{step_status=FAILEDCG;}}
        else if(step_quality<0){step_status=NEGRATIO;}
        else{step_status=CONTRACT;}}

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
            radius=norm_sk_scaled*contract_factor;
            //radius*=contract_factor;
            break;
        case EXPAND:
            radius=std::max(expand_factor*norm_sk_scaled,radius);
            //radius*=expand_factor;
            break;
    };
    return step_status;
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
    p_norm_gk=preconditioner._norm(gk,wd);

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
            reason<<"Negative curvature";
            break;}

        aj=rj.dot(yj)/dBd;
        zj_old=zj;
        zj.noalias()+=aj*dj;

        if(preconditioner._norm(zj,wd)>=radius){
            // find tau>=0 s.t. p intersects trust region
            tau=Find_Tau(zj_old,dj);
            pk.noalias()=zj_old+tau*dj;
            num_CG_iterations=j+1;
            reason<<"Intersect TR bound";
            break;}

        dot_ry=rj.dot(yj);
        rj.noalias()-=aj*(hessian.template selfadjointView<Lower>()*dj).eval();
        
        if(preconditioner._norm(rj,wd)/p_norm_gk<tol){
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
    X._multiply(v,out);
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
GENERIC_TYPE_DEFINITION(TRUST_REGION)
DEFINE_AND_REGISTER_PARSER(TRUST_REGION,void)
{
    auto step=std::make_shared<TRUST_REGION<TV>>();
    step->equation=new NONLINEAR_EQUATION<TV>();
    Parse_String(node["name"],step->name);
    Parse_Scalar(node["precision"],step->precision,step->precision);
    simulation.evolution.push_back(step);
    return 0;
}
