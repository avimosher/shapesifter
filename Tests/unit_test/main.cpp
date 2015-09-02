#define CATCH_CONFIG_MAIN
#include <Data/RIGID_STRUCTURE.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Utilities/RANDOM.h>
#include <catch.hpp>

using namespace Mechanics;
typedef double T;
typedef Eigen::Matrix<T,3,1> TV;

unsigned int Factorial(int number){
    return number<=1?number:Factorial(number-1)*number;
}


TEST_CASE("Factorials are computed","[factorial]"){
    REQUIRE(Factorial(1)==1);
    REQUIRE(Factorial(2)==2);
    REQUIRE(Factorial(3)==6);
    REQUIRE(Factorial(10)==3628800);
}

TEST_CASE("Rigid structure","[rigid structure]"){
    /*auto structure=std::make_shared<RIGID_STRUCTURE<TV>>();
    SECTION("test stepping"){
    }*/

}

/*Estimate_Derivative()
{
    for(int i=1;i<=NTAB;i++){
        hh/=CON;
        a[0][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
        fac=sqr(CON);
        for(int j=1;j<=NTAB;j++){
            a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
            fac=sqr(CON)*fac;
            errt=std::max(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
            if(errt<=*err){
                *err=errt;
                ans=a[j][i];
            }
            if(std::abs(a[i][i]-a[i-1][i-1])>=SAFE*(*err)){break;}
        }
    }
    }*/

TEST_CASE("Derivatives","[derivatives]"){
    RANDOM<T> random;
    int tests=10;
    SECTION("dForce_dVelocity"){
        for(int i=0;i<tests;i++){
            TV x1=random.template Direction<TV>();
            TV x2=random.template Direction<TV>();
            TV relative_position=x2-x1;
            TV direction=relative_position.normalized();
            Eigen::Matrix<T,3,3> derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dForce_dVelocity(relative_position);
            T epsilon=1e-6;
            TV delta=epsilon*random.template Direction<TV>();
            TV estimated_direction=direction+derivative*delta;
            TV actual_direction=(x2+delta-x1).normalized();
            std::cout<<"Quality: "<<(actual_direction-estimated_direction).norm()/epsilon<<std::endl;
            REQUIRE((actual_direction-estimated_direction).norm()<epsilon);
        }
    }
    SECTION("Rotated offset"){
        T epsilon=1e-6;
        for(int i=0;i<tests;i++){
            TV base_offset=random.template Direction<TV>();
            TV spin=random.template Direction<TV>();
            
        }
    }
    SECTION("dRotatedOffset_dSpin"){
        T epsilon=1e-6;
        for(int i=0;i<tests;i++){
            TV base_offset=random.template Direction<TV>();
            TV spin=random.template Direction<TV>();
            Eigen::Matrix<T,3,3> derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spin,ROTATION<TV>::From_Rotation_Vector(spin)*base_offset);
            TV delta=epsilon*random.template Direction<TV>();
            TV estimated_offset=ROTATION<TV>::From_Rotation_Vector(spin)*base_offset+derivative*delta;
            TV actual_offset=ROTATION<TV>::From_Rotation_Vector(spin+delta)*base_offset;
            std::cout<<"Step size: "<<(derivative*delta).norm()<<std::endl;
            std::cout<<"Quality: "<<(actual_offset-estimated_offset).norm()/delta.norm()<<std::endl;
            REQUIRE((actual_offset-estimated_offset).norm()<epsilon/2);
        }
    }
    SECTION("dForce_dSpin"){
        for(int i=0;i<tests;i++){
            std::vector<TV> positions(2);
            std::vector<TV> base_offsets(2);
            std::vector<TV> spins(2);
            std::vector<ROTATION<TV>> rotations(2);
            std::cout<<"TEST"<<std::endl;
            for(int j=0;j<2;j++){
                positions[j]=random.template Direction<TV>();
                base_offsets[j]=random.template Direction<TV>();
                spins[j]=random.template Direction<TV>();
                rotations[j]=ROTATION<TV>::From_Rotation_Vector(spins[j]);
                std::cout<<"Position: "<<positions[j].transpose()<<std::endl;
                std::cout<<"Base offset: "<<base_offsets[j].transpose()<<std::endl;
                std::cout<<"Spins: "<<spins[j].transpose()<<std::endl;
            }
            TV relative_position=positions[1]+rotations[1]*base_offsets[1]-(positions[0]+rotations[0]*base_offsets[0]);
            TV direction=relative_position.normalized();
            T epsilon=1e-8;
            for(int s1=0;s1<2;s1++){
                int overall_sign=s1==0?1:-1;
                for(int s2=0;s2<2;s2++){
                    int term_sign=s1==s2?1:-1;
                    Eigen::Matrix<T,3,3> derivative=term_sign*RIGID_STRUCTURE_INDEX_MAP<TV>::dForce_dSpin(relative_position,spins[s2],rotations[s2]*base_offsets[s2]);
                    //Eigen::Matrix<T,1,3> derivative=relative_position.transpose()*RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spins[s2],rotations[s2]*base_offsets[s2]);
                    //T estimated_norm=norm+derivative*delta;
                    TV delta=epsilon*random.template Direction<TV>();
                    TV estimated_direction=overall_sign*(direction)+derivative*delta;
                    std::vector<TV> mod_spins(2);
                    for(int j=0;j<2;j++){mod_spins[j]=spins[j];}
                    mod_spins[s2]+=delta;
                    TV actual_direction=overall_sign*(positions[1]+ROTATION<TV>::From_Rotation_Vector(mod_spins[1])*base_offsets[1]-(positions[0]+ROTATION<TV>::From_Rotation_Vector(mod_spins[0])*base_offsets[0])).normalized();
                    std::cout<<"Direction: "<<direction.transpose()<<" estimated: "<<estimated_direction.transpose()<<" actual: "<<actual_direction.transpose()<<" relative: "<<relative_position.transpose()<<std::endl;
                    std::cout<<"Quality: "<<(actual_direction-estimated_direction).norm()/delta.norm()<<std::endl;
                    REQUIRE((actual_direction-estimated_direction).norm()<2*delta.norm());
                }
            }
        }
    }
}
