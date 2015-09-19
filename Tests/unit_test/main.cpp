#define CATCH_CONFIG_MAIN
#include <Data/RIGID_STRUCTURE.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Equation/NONLINEAR_EQUATION.h>
#include <Force/FORCE.h>
#include <Force/VOLUME_EXCLUSION_CONSTRAINT.h>
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

TEST_CASE("VOLUME_EXCLUSION_CONSTRAINT","[derivatives]"){
    RANDOM<T> random;
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
    }
}

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
            //std::cout<<"Quality: "<<(actual_direction-estimated_direction).norm()/epsilon<<std::endl;
            REQUIRE((actual_direction-estimated_direction).norm()<epsilon);
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
            //std::cout<<"Step size: "<<(derivative*delta).norm()<<std::endl;
            //std::cout<<"Quality: "<<(actual_offset-estimated_offset).norm()/delta.norm()<<std::endl;
            REQUIRE((actual_offset-estimated_offset).norm()<epsilon/2);
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
                    Eigen::Matrix<T,3,3> derivative=term_sign*RIGID_STRUCTURE_INDEX_MAP<TV>::dForce_dSpin(relative_position,spins[s2],rotations[s2]*base_offsets[s2]);
                    //Eigen::Matrix<T,1,3> derivative=relative_position.transpose()*RIGID_STRUCTURE_INDEX_MAP<TV>::dRotatedOffset_dSpin(spins[s2],rotations[s2]*base_offsets[s2]);
                    //T estimated_norm=norm+derivative*delta;
                    TV delta=epsilon*random.template Direction<TV>();
                    TV estimated_direction=overall_sign*(direction)+derivative*delta;
                    std::vector<TV> mod_spins(2);
                    for(int j=0;j<2;j++){mod_spins[j]=spins[j];}
                    mod_spins[s2]+=delta;
                    TV actual_direction=overall_sign*(positions[1]+ROTATION<TV>::From_Rotation_Vector(mod_spins[1])*base_offsets[1]-(positions[0]+ROTATION<TV>::From_Rotation_Vector(mod_spins[0])*base_offsets[0])).normalized();
                    //std::cout<<"Direction: "<<direction.transpose()<<" estimated: "<<estimated_direction.transpose()<<" actual: "<<actual_direction.transpose()<<" relative: "<<relative_position.transpose()<<std::endl;
                    std::cout<<"Quality: "<<(actual_direction-estimated_direction).norm()/delta.norm()<<std::endl;
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
            Eigen::Matrix<T,3,3> derivative=RIGID_STRUCTURE_INDEX_MAP<TV>::dPenaltyForce_dVelocity(relative_position,threshold);
            TV delta=epsilon*random.template Direction<TV>();
            TV estimated_force=force+derivative*delta;
            TV new_relative=(positions[1]+delta-positions[0]);
            TV actual_force=sqr(new_relative.norm()-threshold)*new_relative.normalized();
            std::cout<<"Quality: "<<(actual_force-estimated_force).norm()/delta.norm()<<std::endl;
            REQUIRE((actual_force-estimated_force).norm()<delta.norm());
        }
    }
}
