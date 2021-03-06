#define CATCH_CONFIG_MAIN
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Force/FORCE.h>
#include <Force/VOLUME_EXCLUSION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Math/F.h>
#include <Math/F_NF.h>
#include <Math/NF.h>
#include <Math/NFINV.h>
#include <Math/R1XRCXR2INV.h>
#include <Math/RCF_NF.h>
#include <Math/RXO.h>
#include <Math/Relative_Position_Force.h>
#include <Math/Spring_Force.h>
#include <Utilities/RANDOM.h>
#include <Eigen/CXX11/Tensor>
#include <Eigen/KroneckerProduct>
#include <catch.hpp>

using namespace Mechanics;
using namespace Eigen;
typedef double T;
typedef Matrix<T,3,1> TV;
typedef Matrix<T,3,1> T_SPIN;
typedef Matrix<T,3,3> M_VxV;
typedef TensorFixedSize<T,Sizes<3,3,3>> T_TENSOR;
typedef Dimension::LINEARITY LINEARITY;

M_VxV Contract(const T_TENSOR& t,const TV& v,const std::array<int,3>& indices){
    M_VxV result;result.setZero();
    std::array<int,3> index{};
    for(index[0]=0;index[0]<3;index[0]++){
        for(index[1]=0;index[1]<3;index[1]++){
            for(index[2]=0;index[2]<3;index[2]++){
                result(index[1],index[2])+=t(index[indices[0]],index[indices[1]],index[indices[2]])*v(index[0]);
            }}}
    return result;
}

TV Contract(const T_TENSOR& t,const TV& v1,const TV& v2,const std::array<int,3>& indices){
    TV result;result.setZero();
    std::array<int,3> index{};
    for(index[0]=0;index[0]<3;index[0]++){
        for(index[1]=0;index[1]<3;index[1]++){
            for(index[2]=0;index[2]<3;index[2]++){
                result(index[2])+=t(index[indices[0]],index[indices[1]],index[indices[2]])*v1(index[0])*v2(index[1]);
            }}}
    return result;
}

TV Evaluate(const std::array<TV,2>& positions,const std::array<T_SPIN,2>& spins,const std::array<TV,2>& offsets){
    return positions[1]+ROTATION<TV>::From_Rotation_Vector(spins[1])*offsets[1]-
        (positions[0]+ROTATION<TV>::From_Rotation_Vector(spins[0])*offsets[0]);}

TEST_CASE("Hessian"){
    RANDOM<T> random;
    T epsilon=1e-3;
    T divisor=2;
    std::array<TV,2> positions;
    std::array<T_SPIN,2> spins;
    std::array<TV,2> offsets,spun_offsets;
    for(int i=0;i<2;i++){
        positions[i]=random.template Direction<TV>();
        spins[i]=random.template Direction<T_SPIN>();
        offsets[i]=random.template Direction<TV>();
        spun_offsets[i]=ROTATION<TV>::From_Rotation_Vector(spins[i])*offsets[i];}
    TV f=Evaluate(positions,spins,offsets);
    TV dx=random.template Direction<TV>();

    SECTION("dnf_dVelocity"){
        TV dnf_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::dnf_dVelocity<LINEARITY::LINEAR>(f,f.norm(),1,spins[0],offsets[0]);
        auto testlambda=[&](T eps){
            T predicted=dnf_dv.dot(eps*dx);
            T actual=Evaluate({positions[0],positions[1]+eps*dx},spins,offsets).norm()-f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);}

    SECTION("dnf_dSpin"){
        T_SPIN spin=random.template Direction<T_SPIN>();
        TV dnf_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::dnf_dVelocity<LINEARITY::ANGULAR>(f,f.norm(),1,spins[1],offsets[1]);
        auto testlambda=[&](T eps){
            T predicted=dnf_ds.dot(eps*dx);
            T actual=Evaluate(positions,{spins[0],spins[1]+eps*dx},offsets).norm()-f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);

        M_VxV d2nf_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2n_dVelocity2<LINEARITY::ANGULAR,LINEARITY::ANGULAR>(f,{1,1},{spins[1],spins[1]},{offsets[1],offsets[1]});
        auto test_second=[&](T eps){
            T predicted=dnf_ds.dot(eps*dx)+(T).5*eps*eps*dx.transpose()*d2nf_ds2*dx;
            T actual=Evaluate(positions,{spins[0],spins[1]+eps*dx},offsets).norm()-f.norm();
            return actual-predicted;};
        ratio=test_second(epsilon)/test_second(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("dnfinv_dVelocity","d2nfinv_dVelocity2"){
        TV dnfinv_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::dnfinv_dVelocity<LINEARITY::LINEAR>(f,f.norm(),1,spins[1],offsets[1]);
        M_VxV d2nfinv_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::d2nfinv_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(f,f.norm(),{1,1},spins,offsets);
        auto testlambda=[&](T eps){
            T predicted=dnfinv_dv.dot(eps*dx)+((T).5*eps*eps*dx.transpose()*d2nfinv_dv*dx);
            T actual=1/Evaluate({positions[0],positions[1]+eps*dx},spins,offsets).norm()-1/f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("dnfinv_dSpin"){
        TV dnfinv_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::dnfinv_dVelocity<LINEARITY::ANGULAR>(f,f.norm(),1,spins[1],offsets[1]);
        auto testlambda=[&](T eps){
            T predicted=dnfinv_dv.dot(eps*dx);
            T actual=1/Evaluate(positions,{spins[0],spins[1]+eps*dx},offsets).norm()-1/f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);}

    SECTION("d2nfinv_dSpin2"){
        M_VxV d2nfinv_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::d2nfinv_dVelocity2<LINEARITY::ANGULAR,LINEARITY::ANGULAR>(f,f.norm(),{1,1},{spins[1],spins[1]},{offsets[1],offsets[1]});
        TV d0=RIGID_STRUCTURE_INDEX_MAP<TV>::dnfinv_dVelocity<LINEARITY::ANGULAR>(f,f.norm(),1,spins[1],offsets[1]);
        auto testlambda=[&](T eps){
            T predicted=d0.dot(eps*dx)+((T).5*eps*eps*dx.transpose()*d2nfinv_dv*dx);
            T actual=1/Evaluate(positions,{spins[0],spins[1]+eps*dx},offsets).norm()-1/f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("dnainv_dA"){
        TV dnainv_da=RIGID_STRUCTURE_INDEX_MAP<TV>::dnainv_dA(f,f.norm());
        auto testlambda=[&](T eps){
            T predicted=dnainv_da.dot(dx)*eps;
            T actual=1/(f+eps*dx).norm()-1/f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);}

    SECTION("da_na_dA"){
        M_VxV da_na_da=RIGID_STRUCTURE_INDEX_MAP<TV>::da_na_dA(f,f.norm());
        T_TENSOR d2a_na_da2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2a_na_dA2(f,f.norm());
        auto test_second=[&](T eps){
            TV contracted=Contract(d2a_na_da2,dx,dx,{0,1,2});
            TV predicted=da_na_da*dx*eps+((T).5*eps*eps*contracted);
            TV actual=(f+eps*dx).normalized()-f.normalized();
            return (actual-predicted).norm();};
        T ratio=test_second(epsilon)/test_second(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("df_nf_dVelocity"){
        M_VxV df_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::df_nf_dVelocity<LINEARITY::LINEAR>(f,1,spins[1],offsets[1]);
        auto testlambda=[&](T eps){
            TV predicted=df_dv*(eps*dx);
            TV actual=Evaluate({positions[0],positions[1]+eps*dx},spins,offsets).normalized()-f.normalized();
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);}

    SECTION("df_dSpin"){
        M_VxV df_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::df_dVelocity<LINEARITY::ANGULAR>(1,spins[1],offsets[1]);
        auto testlambda=[&](T eps){
            TV predicted=df_ds*(eps*dx);
            TV actual=Evaluate(positions,{spins[0],spins[1]+eps*dx},offsets)-Evaluate(positions,spins,offsets);
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);}

    SECTION("d2f_nf_dVelocity2"){
        T_TENSOR d2f_dv2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2f_nf_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(f,{1,1},spins,offsets);
        M_VxV df_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::df_nf_dVelocity<LINEARITY::LINEAR>(f,1,spins[1],offsets[1]);
        auto testlambda=[&](T eps){
            TV predicted=df_dv*(eps*dx)+(T).5*eps*eps*Contract(d2f_dv2,dx,dx,{0,1,2});
            TV actual=Evaluate({positions[0],positions[1]+eps*dx},spins,offsets).normalized()-f.normalized();
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("d2f_nf_dSpin2"){
        T_TENSOR d2f_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2f_nf_dVelocity2<LINEARITY::ANGULAR,LINEARITY::ANGULAR>(f,{1,1},{spins[1],spins[1]},{offsets[1],offsets[1]});
        M_VxV df_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::df_nf_dVelocity<LINEARITY::ANGULAR>(f,1,spins[1],offsets[1]);
        auto testlambda=[&](T eps){
            TV predicted=df_ds*(eps*dx)+(T).5*eps*eps*Contract(d2f_ds2,dx,dx,{0,1,2});
            TV actual=Evaluate(positions,{spins[0],spins[1]+eps*dx},offsets).normalized()-f.normalized();
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("d2n_dVelocity2"){
        M_VxV d2n_dv2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2n_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(f,{1,1},spins,offsets);
        TV dnf_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::dnf_dVelocity<LINEARITY::LINEAR>(f,f.norm(),1,spins[1],offsets[1]);

        auto testlambda=[&](T eps){
            T predicted=dnf_dv.dot(eps*dx)+(T)0.5*eps*eps*dx.transpose()*d2n_dv2*dx;
            T actual=Evaluate({positions[0],positions[1]+eps*dx},spins,offsets).norm()-f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("d2n_dSpin2"){
        M_VxV d2n_dv2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2n_dVelocity2<LINEARITY::ANGULAR,LINEARITY::ANGULAR>(f,{1,1},{spins[1],spins[1]},{offsets[1],offsets[1]});
        TV dnf_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::dnf_dVelocity<LINEARITY::ANGULAR>(f,f.norm(),1,spins[1],offsets[1]);

        auto testlambda=[&](T eps){
            T predicted=dnf_dv.dot(eps*dx)+(T)0.5*eps*eps*dx.transpose()*d2n_dv2*dx;
            T actual=Evaluate(positions,{spins[0],spins[1]+eps*dx},offsets).norm()-f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    /*
      for a proper evaluation, in theory, have to vary all of the things: which derivatives wrt, 

      for a scalar function:
      assume we have 1st and second derivatives for it, wrt angular and linear
      build all of the necessary first derivatives (four) and second derivatives (16)
      


     */

    // for vector functions

    // VTYPE: LINEAR or ANGULAR
    // VSIGN:: -1 or 1
    // V: 0 or 1 (storage index)
    //
    // need second derivatives for all combinations (4 x 4)

    std::array<std::array<TV,2>,2> dxs;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            dxs[i][j]=random.template Direction<TV>();
        }}

    SECTION("F"){
        T ratio=F<TV>::Test_Error(positions,spins,offsets,dxs,epsilon)/F<TV>::Test_Error(positions,spins,offsets,dxs,epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("NF"){
        T ratio=NF<TV>::Test_Error(positions,spins,offsets,dxs,epsilon)/NF<TV>::Test_Error(positions,spins,offsets,dxs,epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("NFINV"){
        T ratio=NFINV<TV>::Test_Error(positions,spins,offsets,dxs,epsilon)/NFINV<TV>::Test_Error(positions,spins,offsets,dxs,epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("F_NF"){
        T ratio=F_NF<TV>::Test_Error(positions,spins,offsets,dxs,epsilon)/F_NF<TV>::Test_Error(positions,spins,offsets,dxs,epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("RXO"){
        T ratio=RXO<TV,1>::Test_Error(positions,spins,offsets,dxs,epsilon)/RXO<TV,1>::Test_Error(positions,spins,offsets,dxs,epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    SECTION("RXF_NF"){
        T ratio0=RCF_NF<TV,0>::Test_Error(positions,spins,offsets,dxs,epsilon)/RCF_NF<TV,0>::Test_Error(positions,spins,offsets,dxs,epsilon/divisor);
        REQUIRE(fabs(ratio0-cube(divisor))<0.1);
    }

    SECTION("RXF_NF"){
        T ratio1=RCF_NF<TV,1>::Test_Error(positions,spins,offsets,dxs,epsilon)/RCF_NF<TV,1>::Test_Error(positions,spins,offsets,dxs,epsilon/divisor);
        REQUIRE(fabs(ratio1-cube(divisor))<0.1);
    }

    SECTION("practical"){
        std::vector<Triplet<T>> hessian_terms;
        TV f=F<TV>::Evaluate(positions,spins,offsets);
        Relative_Position_Force<TV>::Second_Derivatives(f,spins,offsets,{0,1},dxs,hessian_terms);
    }

    SECTION("Spring_Force"){
        T target=3;
        T stiffness=10;
        M_VxV derivative=Spring_Force<TV>::template First_Derivative<1,1,0,1>(stiffness,target,f,spins,offsets);
        auto testlambda=[&](T eps){
            TV predicted=derivative.transpose()*dxs[0][0]*eps;
            TV actual=Spring_Force<TV>::template Evaluate<1,1>(stiffness,target,{positions[0],positions[1]},{spins[0]+eps*dxs[0][0],spins[1]},offsets)-Spring_Force<TV>::template Evaluate<1,1>(stiffness,target,positions,spins,offsets);
            return (actual-predicted).norm();
        };
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);
    }

    SECTION("association dissociation"){
        ROTATION<TV> RC=ROTATION<TV>::From_Rotation_Vector(positions[0]);
        M_VxV derivative=R1XRCXR2INV<TV>::template First_Derivative<0>(RC,spins);
        auto testlambda=[&](T eps){
            TV predicted=derivative.transpose()*dxs[0][0]*eps;
            TV actual=R1XRCXR2INV<TV>::Evaluate(RC,{spins[0]+eps*dxs[0][0],spins[1]})-R1XRCXR2INV<TV>::Evaluate(RC,spins);
            return (actual-predicted).norm();
        };
        T ratio=testlambda(epsilon)/testlambda(epsilon/2);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);
    }

    SECTION("d2f_nf_dVelocity2 full"){
        std::array<M_VxV,2> df_dvs;
        std::array<TV,2> dxs;
        Matrix<T_TENSOR,2,2> d2f_dv2s;
        TV zero;zero.setZero();
        TV f=Evaluate(positions,spins,{zero,zero});
        for(int s1=0,s1_sgn=-1;s1<2;s1++,s1_sgn+=2){
            for(int s2=0,s2_sgn=-1;s2<2;s2++,s2_sgn+=2){
                d2f_dv2s(s1,s2)=RIGID_STRUCTURE_INDEX_MAP<TV>::d2f_nf_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(f,{s1_sgn,s2_sgn},{spins[s1],spins[s1]},{offsets[s1],offsets[s2]});}
            df_dvs[s1]=RIGID_STRUCTURE_INDEX_MAP<TV>::df_nf_dVelocity<LINEARITY::LINEAR>(f,s1_sgn,spins[s1],offsets[s1]);
            dxs[s1]=random.template Direction<TV>();}

        auto testlambda=[&](T eps){
            TV predicted;predicted.setZero();
            for(int s1=0;s1<2;s1++){
                predicted+=df_dvs[s1].transpose()*(eps*dxs[s1]);
                for(int s2=0;s2<2;s2++){
                    predicted+=.5*eps*eps*dxs[s1].transpose()*Contract(d2f_dv2s(s1,s2),dxs[s2],{2,0,1});}}
            TV actual=Evaluate({positions[0]+eps*dxs[0],positions[1]+eps*dxs[1]},spins,{zero,zero}).normalized()-f.normalized();
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    // ROTATIONAL PARTS
    SECTION("dw_dSpin","d2w_dSpin2"){
        T_SPIN spin=random.template Direction<T_SPIN>();
        T norm_spin=spin.norm();
        T initial=cos(norm_spin/2);
        T_SPIN dw_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::dw_dSpin(spin,norm_spin);
        T_SPIN ds=random.template Direction<T_SPIN>();
        M_VxV d2d_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2w_dSpin2(spin,norm_spin);
        auto testlambda_hessian=[&](T eps){
            T predicted=dw_ds.dot(eps*ds)+(T).5*eps*eps*ds.transpose()*d2d_ds2*ds;;
            T actual=cos((spin+eps*ds).norm()/2)-initial;
            return actual-predicted;};

        T ratio_hessian=testlambda_hessian(epsilon)/testlambda_hessian(epsilon/divisor);
        REQUIRE(fabs(ratio_hessian-cube(divisor))<0.1);}

    SECTION("dq_dSpin","d2q_dSpin2"){
        T_SPIN spin=random.template Direction<T_SPIN>();
        T norm_spin=spin.norm();
        
        T_SPIN initial=sinc(norm_spin/2)*spin/2;
        M_VxV dq_dspin=RIGID_STRUCTURE_INDEX_MAP<TV>::dq_dSpin(spin/norm_spin,norm_spin);
        T_SPIN ds=random.template Direction<T_SPIN>();

        T_TENSOR d2d_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2q_dSpin2(spin,norm_spin);
        auto testlambda_hessian=[&](T eps){
            T_SPIN predicted=dq_dspin*(eps*ds)+(T).5*eps*eps*Contract(d2d_ds2,ds,ds,{0,1,2});
            T_SPIN final_spin=spin+eps*ds;
            T_SPIN actual=sinc(final_spin.norm()/2)*final_spin/2-initial;
            return (actual-predicted).norm();};

        T ratio_hessian=testlambda_hessian(epsilon)/testlambda_hessian(epsilon/divisor);
        REQUIRE(fabs(ratio_hessian-cube(divisor))<0.1);}

    SECTION("d2sxo_dSpin2"){
        M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spins[1],offsets[1]);
        TV delta=random.template Direction<TV>();
        T_TENSOR d2so_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2so_dSpin2(spins[1],offsets[1]);
        auto test_second=[&](T eps){
            T_SPIN dspin=eps*dx;
            TV predicted=derivative*dspin+(T).5*Contract(d2so_ds2,dspin,dspin,{0,1,2});
            TV final=ROTATION<TV>::From_Rotation_Vector(spins[1]+dspin)*offsets[1];
            return (final-spun_offsets[1]-predicted).norm();};
        T ratio=test_second(epsilon)/test_second(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}

    // r x f, where f=(f_2-f_1)/|f_2-f_1|
    SECTION("dtau_dSpin","d2tau_dSpin2"){
        TV tau_initial=spun_offsets[1].cross(f.normalized());
        M_VxV dtau_da=RIGID_STRUCTURE_INDEX_MAP<TV>::dtau_dA<LINEARITY::ANGULAR>(f,1,spins[1],offsets[1]);
        T_TENSOR d2tau_da2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2tau_dV2<LINEARITY::ANGULAR,LINEARITY::ANGULAR>(f,spun_offsets[1],{1,1},{spins[1],spins[1]},{offsets[1],offsets[1]});
        auto test_second=[&](T eps){
            T_SPIN dspin=eps*dx;
            TV predicted=dtau_da*dspin+(T).5*Contract(d2tau_da2,dspin,dspin,{0,1,2});
            TV r_final=ROTATION<TV>::From_Rotation_Vector(spins[1]+dspin)*offsets[1];
            TV tau_final=r_final.cross(Evaluate(positions,{spins[0],spins[1]+dspin},offsets).normalized());
            T error=(tau_final-tau_initial-predicted).norm();
            return error;};
        T ratio=test_second(epsilon)/test_second(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);}
}


TEST_CASE("ASSOCIATION_DISSOCATION_CONSTRAINT"){
    RANDOM<T> random;
    int tests=10;
    T epsilon=1e-4;
    SECTION("dRdS x F"){
        for(int i=0;i<tests;i++){
            /*T_SPIN s1=random.template Direction<T_SPIN>();
            T_SPIN s2=random.template Direction<T_SPIN>();
            ROTATION<TV> r1(ROTATION<TV>::From_Rotation_Vector(s1));
            ROTATION<TV> r2(ROTATION<TV>::From_Rotation_Vector(s2));
            TV x1=random.template Direction<TV>();

            std::cout<<"R2: "<<(r2.inverse()*x1).transpose()<<std::endl;
            std::cout<<"Complicated: "<<((r1.inverse().toRotationMatrix()*(r2.toRotationMatrix()*r1.toRotationMatrix().inverse()).inverse())*x1).transpose()<<std::endl;*/

            T_SPIN spin=random.template Direction<T_SPIN>();
            TV base_offset=random.template Direction<TV>();
            TV f=random.template Direction<TV>();
            ROTATION<TV> rotation(ROTATION<TV>::From_Rotation_Vector(spin));
            T_SPIN dspin=epsilon*random.template Direction<T_SPIN>();

            TV initial=(rotation*base_offset).cross(f);
            TV final=(ROTATION<TV>::From_Rotation_Vector(spin+dspin)*base_offset).cross(f);

            
            M_VxV dFdS;
            
            M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin,base_offset);
            for(int j=0;j<3;j++){
                dFdS.block<3,1>(0,j)=derivative.block<3,1>(0,j).cross(f);
            }
            TV predicted_delta=dFdS*dspin;
            REQUIRE((final-initial-predicted_delta).norm()<1e-2*epsilon);
        }
    }
    SECTION("dCrdS"){
        // derivative of the rotation constraint
        for(int i=0;i<tests;i++){
            T_SPIN spin=random.template Direction<T_SPIN>();
            ROTATION<TV> base_rotation=ROTATION<TV>::From_Rotation_Vector(random.template Direction<T_SPIN>());
            ROTATION<TV> rotation(ROTATION<TV>::From_Rotation_Vector(spin));
            T_SPIN dspin_direction=random.template Direction<T_SPIN>();
            ROTATION<TV> total_rotation=base_rotation*rotation.inverse();
            TV initial=total_rotation.vec();
            M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::Compute_Simple_Orientation_Constraint_Matrix(rotation,total_rotation,-1);
            auto testlambda=[&](T eps){
                T_SPIN dspin=eps*dspin_direction;
                TV final=(base_rotation*ROTATION<TV>::From_Rotation_Vector(spin+dspin).inverse()).vec();
                TV predicted=derivative*dspin;
                return (final-initial-predicted).norm();
            };
            T ratio=testlambda(epsilon)/testlambda(epsilon/2);
            REQUIRE(fabs(ratio-4)<1e-2);
        }
    }
}

TEST_CASE("VOLUME_EXCLUSION_CONSTRAINT","[derivatives]"){
    /*RANDOM<T> random;
    SECTION("derivative"){
        SIMULATION<TV> simulation;
        auto rigid_data=simulation.data.template Find_Or_Create<RIGID_STRUCTURE_DATA<TV>>();
        auto volume_exclusion_constraint=simulation.force.template Find_Or_Create<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
        auto structure=std::make_shared<RIGID_STRUCTURE<TV>>();
        structure->frame.position=random.template Direction<TV>();
        structure->radius=1;
        structure->collision_radius=1;
        structure->Initialize_Inertia(3.5);
        rigid_data->structures.push_back(structure);

        auto structure2=std::make_shared<RIGID_STRUCTURE<TV>>();
        structure2->frame.position=structure->frame.position+TV::UnitX()*1.8;
        structure2->radius=1;
        structure2->collision_radius=1;
        structure2->Initialize_Inertia(3.5);
        rigid_data->structures.push_back(structure2);
        

        NONLINEAR_EQUATION<TV> equation;
        equation.Linearize(simulation.data,simulation.force,1,0,1);
        //std::cout<<equation.jacobian<<std::endl;
        
        REQUIRE(1);
        }*/
}

TEST_CASE("Derivatives","[derivatives]"){
    RANDOM<T> random;
    int tests=10;
    SECTION("dConstraint_dTwist"){
        // the constraint is |x2-x1|
        for(int i=0;i<tests;i++){
            TV x1=random.template Direction<TV>();
            TV x2=random.template Direction<TV>();
            T_SPIN spin=random.template Direction<T_SPIN>();
            TV base_offset=random.template Direction<TV>();;

            TV rotated_offset=ROTATION<TV>::From_Rotation_Vector(spin)*base_offset;
            TV relative_position=x2+rotated_offset-x1;
            Matrix<T,1,6> derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(spin,base_offset,relative_position);

            T epsilon=1e-5;
            Matrix<T,6,1> dx_direction;random.Direction(dx_direction);

            T initial=(x2+rotated_offset-x1).norm();

            auto testlambda=[&](T eps){
                Matrix<T,6,1> dx=eps*dx_direction;
                T final=(x2+dx.block<3,1>(0,0)+ROTATION<TV>::From_Rotation_Vector(spin+dx.block<3,1>(3,0))*base_offset-x1).norm();
                T predicted=derivative*dx;
                //std::cout<<"final: "<<final<<" initial: "<<initial<<" predicted: "<<predicted<<std::endl;
                return (final-initial-predicted);
            };
            T ratio=testlambda(epsilon)/testlambda(epsilon/2);
            //std::cout<<"dC_dT ratio: "<<ratio<<std::endl;
            REQUIRE(fabs(ratio-4)<0.1);
        }
    }
    SECTION("dForce_dVelocity"){
        for(int i=0;i<tests;i++){
            TV x1=random.template Direction<TV>();
            TV x2=random.template Direction<TV>();
            TV relative_position=x2-x1;
            TV direction=relative_position.normalized();
            M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dForce_dVelocity(relative_position);
            T epsilon=1e-6;
            TV delta=epsilon*random.template Direction<TV>();
            TV estimated_direction=direction+derivative*delta;
            TV actual_direction=(x2+delta-x1).normalized();
            //std::cout<<"Quality: "<<(actual_direction-estimated_direction).norm()/epsilon<<std::endl;
            REQUIRE((actual_direction-estimated_direction).norm()<epsilon);
        }
    }
    SECTION("dRotatedOffset_dSpin"){
        T epsilon=1e-6;
        for(int i=0;i<tests;i++){
            TV base_offset=random.template Direction<TV>();
            TV spin=random.template Direction<TV>();
            TV rotated_offset=ROTATION<TV>::From_Rotation_Vector(spin)*base_offset;
            M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin,base_offset);
            TV delta=random.template Direction<TV>();

            auto testlambda=[&](T eps){
                T_SPIN dspin=eps*delta;
                TV predicted=derivative*dspin;
                TV final=ROTATION<TV>::From_Rotation_Vector(spin+dspin)*base_offset;
                T error=(final-rotated_offset-predicted).norm();
                return error;
            };
            T divisor=2;
            T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
            //std::cout<<"Ratio: "<<ratio<<std::endl;
            REQUIRE(fabs(ratio-sqr(divisor))<0.1);
        }
    }
    SECTION("dForce_dSpin"){
        for(int i=0;i<tests;i++){
            std::vector<TV> positions(2);
            std::vector<TV> base_offsets(2);
            std::vector<TV> spins(2);
            std::vector<ROTATION<TV>> rotations(2);
            for(int j=0;j<2;j++){
                positions[j]=random.template Direction<TV>();
                base_offsets[j]=random.template Direction<TV>();
                spins[j]=random.template Direction<TV>();
                rotations[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]);
            }
            TV relative_position=positions[1]+rotations[1]*base_offsets[1]-(positions[0]+rotations[0]*base_offsets[0]);
            TV direction=relative_position.normalized();
            T epsilon=1e-8;
            for(int s1=0;s1<2;s1++){
                int overall_sign=s1==0?-1:1;
                for(int s2=0;s2<2;s2++){
                    int term_sign=s1==s2?1:-1;
                    M_VxV derivative=term_sign*RIGID_STRUCTURE_INDEX_MAP<TV>::dForce_dSpin(relative_position,spins[s2],base_offsets[s2]);
                    TV delta=epsilon*random.template Direction<TV>();
                    TV estimated_direction=overall_sign*(direction)+derivative*delta;
                    std::vector<TV> mod_spins(2);
                    for(int j=0;j<2;j++){mod_spins[j]=spins[j];}
                    mod_spins[s2]+=delta;
                    TV actual_direction=overall_sign*(positions[1]+ROTATION<TV>::From_Rotation_Vector(mod_spins[1])*base_offsets[1]-(positions[0]+ROTATION<TV>::From_Rotation_Vector(mod_spins[0])*base_offsets[0])).normalized();
                    REQUIRE((actual_direction-estimated_direction).norm()<2*delta.norm());}}}}

    SECTION("Penalty force","dPenaltyForce_dVelocity"){
        for(int i=0;i<tests;i++){
            std::vector<TV> positions(2);
            std::vector<TV> base_offsets(2);
            std::vector<TV> spins(2);
            std::vector<ROTATION<TV>> rotations(2);
            for(int j=0;j<2;j++){
                positions[j]=random.template Direction<TV>();
                base_offsets[j]=random.template Direction<TV>();
                spins[j]=random.template Direction<TV>();
                rotations[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]);}
            T threshold=(positions[1]-positions[0]).norm()+random.Uniform((T)0,(T).1);
            //TV relative_position=positions[1]+rotations[1]*base_offsets[1]-(positions[0]+rotations[0]*base_offsets[0]);
            TV relative_position=positions[1]-positions[0];
            TV direction=relative_position.normalized();
            T epsilon=1e-8;
            TV force=sqr(relative_position.norm()-threshold)*relative_position.normalized();
            M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dPenaltyForce_dVelocity(relative_position,threshold);
            TV delta=epsilon*random.template Direction<TV>();
            TV estimated_force=force+derivative*delta;
            TV new_relative=(positions[1]+delta-positions[0]);
            TV actual_force=sqr(new_relative.norm()-threshold)*new_relative.normalized();
            //std::cout<<"Quality: "<<(actual_force-estimated_force).norm()/delta.norm()<<std::endl;
            REQUIRE((actual_force-estimated_force).norm()<delta.norm());
        }
    }

    SECTION("Penalty force","dPenaltyTorque_dSpin"){
        for(int i=0;i<tests;i++){
            std::vector<TV> positions(2);
            std::vector<TV> base_offsets(2);
            std::vector<TV> spins(2);
            std::vector<ROTATION<TV>> rotations(2);
            for(int j=0;j<2;j++){
                positions[j]=random.template Direction<TV>();
                base_offsets[j]=random.template Direction<TV>();
                spins[j]=random.template Direction<TV>();
                rotations[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]);
            }
            TV relative_position=positions[1]+rotations[1]*base_offsets[1]-(positions[0]+rotations[0]*base_offsets[0]);
            T threshold=relative_position.norm()+random.Uniform((T)0,(T).1);
            std::vector<TV> rotated_offsets={rotations[0]*base_offsets[0],rotations[1]*base_offsets[1]};
            T epsilon=1e-8;
            for(int s1=0;s1<2;s1++){
                int overall_sign=s1==0?-1:1;
                for(int s2=0;s2<2;s2++){
                    int term_sign=s1==s2?1:-1;
                    M_VxV derivative=term_sign*RIGID_STRUCTURE_INDEX_MAP<TV>::dPenaltyTorque_dSpin(relative_position,s1,s2,spins[s2],base_offsets[s1],base_offsets[s2],threshold);
                    TV delta=epsilon*random.template Direction<TV>();

                    TV torque=overall_sign*(ROTATION<TV>::From_Rotation_Vector(spins[s1])*base_offsets[s1]).cross(sqr(relative_position.norm()-threshold)*relative_position.normalized());

                    TV estimated_torque=torque+derivative*delta;
                    std::vector<TV> mod_spins(2);
                    for(int j=0;j<2;j++){mod_spins[j]=spins[j];}
                    mod_spins[s2]+=delta;
                    TV new_relative_position=positions[1]+ROTATION<TV>::From_Rotation_Vector(mod_spins[1])*base_offsets[1]-(positions[0]+ROTATION<TV>::From_Rotation_Vector(mod_spins[0])*base_offsets[0]);
                    TV actual_torque=overall_sign*(ROTATION<TV>::From_Rotation_Vector(mod_spins[s1])*base_offsets[s1]).cross(sqr(new_relative_position.norm()-threshold)*new_relative_position.normalized());
                    //std::cout<<"Quality: "<<(actual_torque-estimated_torque).norm()/delta.norm()<<std::endl;
                    REQUIRE((actual_torque-estimated_torque).norm()<delta.norm());

                }
            }
        }        
    }
}
