#define CATCH_CONFIG_MAIN
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Force/FORCE.h>
#include <Force/VOLUME_EXCLUSION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
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

unsigned int Factorial(int number){
    return number<=1?number:Factorial(number-1)*number;
}

M_VxV Contract(const T_TENSOR& t,const Matrix<T,3,1>& v,const std::array<int,3>& indices){
    M_VxV result;result.setZero();
    std::array<int,3> index{};
    for(index[0]=0;index[0]<3;index[0]++){
        for(index[1]=0;index[1]<3;index[1]++){
            for(index[2]=0;index[2]<3;index[2]++){
                result(index[1],index[2])+=t(index[indices[0]],index[indices[1]],index[indices[2]])*v(index[0]);
            }}}
    return result;
}

Matrix<T,3,1> Contract(const T_TENSOR& t,const Matrix<T,3,1>& v1,const Matrix<T,3,1>& v2,const std::array<int,3>& indices){
    Matrix<T,3,1> result;result.setZero();
    std::array<int,3> index{};
    for(index[0]=0;index[0]<3;index[0]++){
        for(index[1]=0;index[1]<3;index[1]++){
            for(index[2]=0;index[2]<3;index[2]++){
                result(index[2])+=t(index[indices[0]],index[indices[1]],index[indices[2]])*v1(index[0])*v2(index[1]);
            }}}
    return result;
}

T Evaluate(const TV& v){
    return v.normalized().sum();
}

TV Evaluate_Vector(const TV& v){
    return v.normalized();
}

TEST_CASE("Hessian"){
    RANDOM<T> random;
    T epsilon=1e-3;
    TV x1=random.template Direction<TV>();
    TV x2=random.template Direction<TV>();
    TV f0=x2-x1;
    TV dx1=random.template Direction<TV>();
    TV dx2=random.template Direction<TV>();
    T divisor=2;
    std::array<T_SPIN,2> spins;
    std::array<TV,2> offsets;

    SECTION("dnf_dVelocity"){
        TV dnf_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::dnf_dVelocity<LINEARITY::LINEAR>(f0,f0.norm(),1,spins[0],offsets[0]);
        auto testlambda=[&](T eps){
            T predicted=dnf_dv.dot(eps*dx2);
            T actual=(x2+eps*dx2-x1).norm()-f0.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);
    }

    SECTION("dnf_dSpin"){
        T_SPIN spin=random.template Direction<T_SPIN>();
        TV offset=random.template Direction<TV>();
        TV spun_offset=ROTATION<TV>::From_Rotation_Vector(spin)*offset;
        TV f=x2+spun_offset-x1;
        TV dnf_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::dnf_dVelocity<LINEARITY::ANGULAR>(f,f.norm(),1,spin,spun_offset);
        auto testlambda=[&](T eps){
            T predicted=dnf_ds.dot(eps*dx2);
            T actual=(x2+ROTATION<TV>::From_Rotation_Vector(spin+eps*dx2)*offset-x1).norm()-f.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);

        M_VxV d2nf_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2n_dVelocity2<LINEARITY::ANGULAR,LINEARITY::ANGULAR>(f,{1,1},{spin,spin},{spun_offset,spun_offset});
        auto test_second=[&](T eps){
            T predicted=dnf_ds.dot(eps*dx2)+(T).5*eps*eps*dx2.transpose()*d2nf_ds2*dx2;
            T actual=(x2+ROTATION<TV>::From_Rotation_Vector(spin+eps*dx2)*offset-x1).norm()-f.norm();
            return actual-predicted;};
        ratio=test_second(epsilon)/test_second(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);
    }

    SECTION("dnfinv_dVelocity"){
        TV dnfinv_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::dnfinv_dVelocity<LINEARITY::LINEAR>(f0,f0.norm(),1,spins[0],offsets[0]);
        auto testlambda=[&](T eps){
            T predicted=dnfinv_dv.dot(eps*dx2);
            T actual=1/(x2+eps*dx2-x1).norm()-1/f0.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);
    }

    SECTION("d2nfinv_dVelocity2"){
        M_VxV d2nfinv_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::d2nfinv_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(f0,f0.norm(),{1,1},spins,offsets);
        // analytical first derivatives
        TV d0=RIGID_STRUCTURE_INDEX_MAP<TV>::dnfinv_dVelocity<LINEARITY::LINEAR>(f0,f0.norm(),1,spins[0],offsets[0]);
        auto testlambda=[&](T eps){
            T predicted=d0.dot(eps*dx2)+((T).5*eps*eps*dx2.transpose()*d2nfinv_dv*dx2);
            T actual=1/(x2+eps*dx2-x1).norm()-1/f0.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);
    }

    SECTION("dnainv_dA"){
        Matrix<T,3,1> dnainv_da=RIGID_STRUCTURE_INDEX_MAP<TV>::dnainv_dA(f0,f0.norm());
        auto testlambda=[&](T eps){
            T predicted=dnainv_da.dot(dx1)*eps;
            T actual=1/(f0+eps*dx1).norm()-1/f0.norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);
    }

    SECTION("da_na_dA"){
        M_VxV da_na_da=RIGID_STRUCTURE_INDEX_MAP<TV>::da_na_dA(f0,f0.norm());
        T_TENSOR d2a_na_da2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2a_na_dA2(f0,f0.norm());
        auto testlambda=[&](T eps){
            TV predicted=da_na_da*dx1*eps;
            TV actual=(f0+eps*dx1).normalized()-f0.normalized();
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);

        auto test_second=[&](T eps){
            TV contracted=Contract(d2a_na_da2,dx1,dx1,{0,1,2});
            TV predicted=da_na_da*dx1*eps+((T).5*eps*eps*contracted);
            TV actual=(f0+eps*dx1).normalized()-f0.normalized();
            return (actual-predicted).norm();};
        ratio=test_second(epsilon)/test_second(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);
    }

    SECTION("df_dVelocity"){
        M_VxV df_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::df_nf_dVelocity<LINEARITY::LINEAR>(f0,1,spins[0],offsets[0]);
        auto testlambda=[&](T eps){
            TV predicted=df_dv.transpose()*(eps*dx2);
            TV actual=Evaluate_Vector(x2+eps*dx2-x1)-Evaluate_Vector(x2-x1);
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);
    }

    SECTION("df_dSpin"){
        T_SPIN spin=random.template Direction<T_SPIN>();
        TV offset=random.template Direction<TV>();
        TV r=ROTATION<TV>::From_Rotation_Vector(spin)*offset;
        M_VxV df_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::df_dVelocity<LINEARITY::ANGULAR>(1,spin,r);
        auto testlambda=[&](T eps){
            TV predicted=df_ds*(eps*dx2);
            TV actual=(x2+ROTATION<TV>::From_Rotation_Vector(spin+eps*dx2)*offset-x1)-(x2+r-x1);
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);
    }

    SECTION("d2f_dVelocity2"){
        T_TENSOR d2f_dv2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2f_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(f0,{1,1},spins,offsets);
        M_VxV df_dv=RIGID_STRUCTURE_INDEX_MAP<TV>::df_nf_dVelocity<LINEARITY::LINEAR>(f0,1,spins[0],offsets[0]);
        M_VxV delta;delta.setZero();
        auto testlambda=[&](T eps){
            TV predicted=df_dv.transpose()*(eps*dx2);
            predicted+=.5*eps*eps*dx2.transpose()*Contract(d2f_dv2,dx2,{2,1,0});
            TV actual=Evaluate_Vector(x2+eps*dx2-x1)-Evaluate_Vector(x2-x1);
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);
    }

    SECTION("d2n_dVelocity2"){
        M_VxV d2n_dv2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2n_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(f0,{1,1},spins,offsets);
        Matrix<T,3,1> d0=RIGID_STRUCTURE_INDEX_MAP<TV>::dnf_dVelocity<LINEARITY::LINEAR>(f0,f0.norm(),1,spins[0],offsets[0]);

        auto testlambda=[&](T eps){
            T predicted=d0.dot(eps*dx2)+(T)0.5*eps*eps*dx2.transpose()*d2n_dv2*dx2;
            T actual=(x2+eps*dx2-x1).norm()-(x2-x1).norm();
            return actual-predicted;};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);
    }

    SECTION("d2n_dVelocity2 full"){
        std::array<M_VxV,2> df_dvs;
        Matrix<T_TENSOR,2,2> d2f_dv2s;
        for(int s1=0,s1_sgn=-1;s1<2;s1++,s1_sgn+=2){
            for(int s2=0,s2_sgn=-1;s2<2;s2++,s2_sgn+=2){
                d2f_dv2s(s1,s2)=RIGID_STRUCTURE_INDEX_MAP<TV>::d2f_dVelocity2<LINEARITY::LINEAR,LINEARITY::LINEAR>(f0,{s1_sgn,s2_sgn},{spins[s1],spins[s1]},{offsets[s1],offsets[s2]});}
            df_dvs[s1]=RIGID_STRUCTURE_INDEX_MAP<TV>::df_nf_dVelocity<LINEARITY::LINEAR>(f0,s1_sgn,spins[s1],offsets[s1]);
        }
        std::array<TV,2> dxs={dx1,dx2};

        auto testlambda=[&](T eps){
            TV predicted;predicted.setZero();
            for(int s1=0,s1_sgn=-1;s1<2;s1++,s1_sgn+=2){
                predicted+=df_dvs[s1].transpose()*(eps*dxs[s1]);
                for(int s2=0,s2_sgn=-1;s2<2;s2++,s2_sgn+=2){
                    predicted+=.5*eps*eps*dxs[s1].transpose()*Contract(d2f_dv2s(s1,s2),dxs[s2],{2,0,1});}}
            TV actual=Evaluate_Vector(x2+eps*dxs[1]-x1-eps*dxs[0])-Evaluate_Vector(x2-x1);
            return (actual-predicted).norm();};
        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-cube(divisor))<0.1);
    }

    // ROTATIONAL PARTS
    SECTION("dw_dSpin","d2w_dSpin2"){
        T_SPIN spin=random.template Direction<T_SPIN>();
        T norm_spin=spin.norm();
        
        T initial=cos(norm_spin/2);
        T_SPIN dw_ds=RIGID_STRUCTURE_INDEX_MAP<TV>::dw_dSpin(spin,norm_spin);
        T_SPIN ds=random.template Direction<T_SPIN>();

        auto testlambda=[&](T eps){
            T predicted=dw_ds.dot(eps*ds);
            T actual=cos((spin+eps*ds).norm()/2)-initial;
            return actual-predicted;};

        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);

        M_VxV d2d_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2w_dSpin2(spin,norm_spin);
        auto testlambda_hessian=[&](T eps){
            T predicted=dw_ds.dot(eps*ds)+(T).5*eps*eps*ds.transpose()*d2d_ds2*ds;;
            T actual=cos((spin+eps*ds).norm()/2)-initial;
            return actual-predicted;};

        T ratio_hessian=testlambda_hessian(epsilon)/testlambda_hessian(epsilon/divisor);
        REQUIRE(fabs(ratio_hessian-cube(divisor))<0.1);
    }

    SECTION("dq_dSpin","d2q_dSpin2"){
        T_SPIN spin=random.template Direction<T_SPIN>();
        T norm_spin=spin.norm();
        
        T_SPIN initial=sinc(norm_spin/2)*spin/2;
        M_VxV dq_dspin=RIGID_STRUCTURE_INDEX_MAP<TV>::dq_dSpin(spin/norm_spin,norm_spin);
        T_SPIN ds=random.template Direction<T_SPIN>();

        auto testlambda=[&](T eps){
            T_SPIN predicted=dq_dspin*(eps*ds);
            T_SPIN final_spin=spin+eps*ds;
            T_SPIN actual=sinc(final_spin.norm()/2)*final_spin/2-initial;
            return (actual-predicted).norm();};

        T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
        REQUIRE(fabs(ratio-sqr(divisor))<0.1);

        T_TENSOR d2d_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2q_dSpin2(spin,norm_spin);
        auto testlambda_hessian=[&](T eps){
            T_SPIN predicted=dq_dspin*(eps*ds)+(T).5*eps*eps*Contract(d2d_ds2,ds,ds,{0,1,2});
            T_SPIN final_spin=spin+eps*ds;
            T_SPIN actual=sinc(final_spin.norm()/2)*final_spin/2-initial;
            return (actual-predicted).norm();};

        T ratio_hessian=testlambda_hessian(epsilon)/testlambda_hessian(epsilon/divisor);
        REQUIRE(fabs(ratio_hessian-cube(divisor))<0.1);
    }


    SECTION("d2sxo_dSpin2"){
        int tests=10;
        for(int i=0;i<tests;i++){
            TV base_offset=random.template Direction<TV>();
            TV spin=random.template Direction<TV>();
            TV rotated_offset=ROTATION<TV>::From_Rotation_Vector(spin)*base_offset;
            M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin,rotated_offset);
            TV delta=random.template Direction<TV>();

            auto testlambda=[&](T eps){
                T_SPIN dspin=eps*delta;
                TV predicted=derivative*dspin;
                TV final=ROTATION<TV>::From_Rotation_Vector(spin+dspin)*base_offset;
                T error=(final-rotated_offset-predicted).norm();
                return error;
            };
            T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
            REQUIRE(fabs(ratio-sqr(divisor))<0.1);

            T_TENSOR d2so_ds2=RIGID_STRUCTURE_INDEX_MAP<TV>::d2so_dSpin2(spin,rotated_offset);

            auto test_second=[&](T eps){
                T_SPIN dspin=eps*delta;
                TV predicted=derivative*dspin+(T).5*Contract(d2so_ds2,dspin,dspin,{0,1,2});
                TV final=ROTATION<TV>::From_Rotation_Vector(spin+dspin)*base_offset;
                return (final-rotated_offset-predicted).norm();
            };
            ratio=test_second(epsilon)/test_second(epsilon/divisor);
            REQUIRE(fabs(ratio-cube(divisor))<0.1);
        }
    }

    // r x f, where f=(f_2-f_1)/|f_2-f_1|, f_2=x_2+r
    SECTION("dtau_dSpin"){
        int tests=10;
        for(int i=0;i<tests;i++){
            TV base_offset=random.template Direction<TV>();
            TV spin=random.template Direction<TV>();
            TV r=ROTATION<TV>::From_Rotation_Vector(spin)*base_offset;
            TV x1=random.template Direction<TV>();
            TV x2=random.template Direction<TV>();
            TV f=(x2+r-x1);
            TV tau_initial=r.cross(f.normalized());
            M_VxV dr_da=RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin,r);
            M_VxV df_da=RIGID_STRUCTURE_INDEX_MAP<TV>::df_nf_dVelocity<LINEARITY::ANGULAR>(f,1,spin,r);
            TV delta=random.template Direction<TV>();

            auto testlambda=[&](T eps){
                T_SPIN dspin=eps*delta;
                TV predicted=(-Cross_Product_Matrix(f.normalized())*dr_da+Cross_Product_Matrix(r)*df_da)*dspin;
                //TV predicted=df_da*dspin;
                TV r_final=ROTATION<TV>::From_Rotation_Vector(spin+dspin)*base_offset;
                TV tau_final=r_final.cross((x2+r_final-x1).normalized());
                T error=(tau_final-tau_initial-predicted).norm();
                //T error=((x2+r_final-x1).normalized()-f.normalized()-predicted).norm();
                return error;
            };
            T ratio=testlambda(epsilon)/testlambda(epsilon/divisor);
            REQUIRE(fabs(ratio-sqr(divisor))<0.1);
        }
    }
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
            
            M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin,ROTATION<TV>::From_Rotation_Vector(spin)*base_offset);
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
            Matrix<T,1,6> derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dConstraint_dTwist(spin,rotated_offset,relative_position);

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
            M_VxV derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin,rotated_offset);
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
                //std::cout<<"Position: "<<positions[j].transpose()<<std::endl;
                //std::cout<<"Base offset: "<<base_offsets[j].transpose()<<std::endl;
                //std::cout<<"Spins: "<<spins[j].transpose()<<std::endl;
            }
            TV relative_position=positions[1]+rotations[1]*base_offsets[1]-(positions[0]+rotations[0]*base_offsets[0]);
            TV direction=relative_position.normalized();
            T epsilon=1e-8;
            for(int s1=0;s1<2;s1++){
                int overall_sign=s1==0?-1:1;
                for(int s2=0;s2<2;s2++){
                    int term_sign=s1==s2?1:-1;
                    M_VxV derivative=term_sign*RIGID_STRUCTURE_INDEX_MAP<TV>::dForce_dSpin(relative_position,spins[s2],rotations[s2]*base_offsets[s2]);
                    //Eigen::Matrix<T,1,3> derivative=relative_position.transpose()*RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spins[s2],rotations[s2]*base_offsets[s2]);
                    //T estimated_norm=norm+derivative*delta;
                    TV delta=epsilon*random.template Direction<TV>();
                    TV estimated_direction=overall_sign*(direction)+derivative*delta;
                    std::vector<TV> mod_spins(2);
                    for(int j=0;j<2;j++){mod_spins[j]=spins[j];}
                    mod_spins[s2]+=delta;
                    TV actual_direction=overall_sign*(positions[1]+ROTATION<TV>::From_Rotation_Vector(mod_spins[1])*base_offsets[1]-(positions[0]+ROTATION<TV>::From_Rotation_Vector(mod_spins[0])*base_offsets[0])).normalized();
                    //std::cout<<"Direction: "<<direction.transpose()<<" estimated: "<<estimated_direction.transpose()<<" actual: "<<actual_direction.transpose()<<" relative: "<<relative_position.transpose()<<std::endl;
                    //std::cout<<"Quality: "<<(actual_direction-estimated_direction).norm()/delta.norm()<<std::endl;
                    REQUIRE((actual_direction-estimated_direction).norm()<2*delta.norm());
                }
            }
        }
    }

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
                rotations[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]);
                //std::cout<<"Position: "<<positions[j].transpose()<<std::endl;
                //std::cout<<"Base offset: "<<base_offsets[j].transpose()<<std::endl;
                //std::cout<<"Spins: "<<spins[j].transpose()<<std::endl;
            }
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
                    M_VxV derivative=term_sign*RIGID_STRUCTURE_INDEX_MAP<TV>::dPenaltyTorque_dSpin(relative_position,s1,s2,spins[s2],rotated_offsets[s1],rotated_offsets[s2],threshold);
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
