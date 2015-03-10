using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
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
Run
{
    // based on trustOptim implementation
    int iteration=0;
    auto status=CONTINUE;
    
    do{
        iteration++;
        status=Update_One_Step();
        if(nrm_gk/sqrt(double(nvars))<=prec){
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
Update_Preconditioner()
{
    bool success=FALSE;
    double alpha,beta,bmin;

    Vector T(nvars);
    SparseMatrix<T> BB(nvars,nvars);
    BB=Bk.template selfadjointView<Lower>();
    
    for(int j=0;i<BB.outerSize();j++){
        T(j)=sqrt(BB.innerVector(j).dot(BB.innerVector(j)));
    }

    for(int j=0;j<BB.outerSize();j++){
        for(typename SparseMatrix<T>::InnerIterator it(BB,j);it;++it){
            BB.coeffRef(it.row(),j)*=1/sqrt(T(it.row())*T(j));
        }
    }

    beta=sqrt(BB.cwiseAbs2().sum());
    bmin=BB.coeff(0,0);
    for(int j=0;j<nvars;j++){
        bmin=std::min(bmin,BBcoeff(j,j));
    }
    
    if(bmin>0){alpha=0;}
    else{alpha=beta/2;}

    int ii=0;
    do{
        ii++;
        PrecondLLt.factorize(BB);
        if(PrecondLLt.info()==Eigen::Success){
            success=TRUE;
        }
        else{
            alpha=std::max(2*alpha,beta/2)-alpha;
            for(int j=0;j<nvars;++){
                BB.coeffRef(j,j)+=alpha;
            }
        }
    }while(!success);
}
///////////////////////////////////////////////////////////////////////
Update_Hessian()
{
    // call to NONLINEAR_EQUATION
    equation->Get_Hessian();
    Bk*=function_scale_factor;
}
///////////////////////////////////////////////////////////////////////
Update_One_Step()
{
    auto step_status=UNKNOWN;
    Solve_Trust_CG(sk);
    norm_sk_scaled=Get_Norm_Sk(PrecondLLt);
    if(!Finite(norm_sk_scaled)){step_status=FAILEDCG;}
    else{
        try_x=xk+sk;
        Get_F(try_x,try_f);
        if(Finite(try_f)){
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
            if(Finite(try_g.norm())){
                try_g*=function_scale_factor;
                yk=try_g-gk;
                f=try_f;
                xk+=sk;
                gk=try_g;
                norm_gk=gk.norm();
                if(ap>expand_threshold_ap && norm_sk_scaled>=expand_threshold_rad*rad){
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
            rad*=contract_factor;
            break;
        case EXPAND:
            rad*=expand_factor;
            break;
    };
    return step_status;
}
///////////////////////////////////////////////////////////////////////
Solve_Trust_CG()
{
    MatrixBase<T>& pk=const_cast<MatrixBase<T>&>(pk_);
    double norm_rj,dot_ry,dot_ry_old,norm_zj,aj,bj,tau,dBd,norm_gk;
    int j;
    double crit;
    
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

        Upz(PrecondLLt,zj,wd);
        norm_zj=wd.norm();

        if(norm_zj>=rad){
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
Finite(T x)
{
    return (std::abs(x)<=__DBL_MAX__) && (x==x);
}
///////////////////////////////////////////////////////////////////////
Get_Norm_Sk(const SimplicialLLT<T>& X)
{
    double res=(X.matrixU()*(X.permutationPinv()*sk).eval()).norm();
    return res;
}
///////////////////////////////////////////////////////////////////////
Get_F()
{
    // TODO: CALL MODEL
    return equation->Evaluate();
}
///////////////////////////////////////////////////////////////////////
Gradient()
{
    // TODO: CALL MODEL
    return equation->Gradient();
}
///////////////////////////////////////////////////////////////////////
UPz(const SimplicialLLT<T>& X,const Vector& v,Vector& out)
{
    out=X.permutationP()*v;
    out=X.matrixU().template triangularView<Upper>()*out;
}
///////////////////////////////////////////////////////////////////////
Find_Tau(const Vector& z,const Vector& d)
{
    UPz(PrecondLLt,d,wd);
    UPz(PrecondLLtz,wz);
    
    double d2=wd.squaredNorm();
    double z2=z.squaredNorm();
    double zd=wd.dot(wz);
    
    double root=zd*zd-d2*(z2-rad*rad);
    double tau=(sqrt(root)-zd)/d2;
    return tau;
}
///////////////////////////////////////////////////////////////////////
Get_FDF()
{
    // TODO: CALL MODEL
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
