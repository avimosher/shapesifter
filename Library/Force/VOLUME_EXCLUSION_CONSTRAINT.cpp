#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Force/FORCE.h>
#include <Force/VOLUME_EXCLUSION_CONSTRAINT.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/HASHING.h>
#include <Utilities/OSG_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <iostream>
#include <osg/Geometry>
#include <math.h>
#include <osg/Geode>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Unpack_Forces(const Matrix<T,Dynamic,1>& forces)
{
    stored_forces=forces;
    for(int i=0;i<constraints.size();i++){
        force_memory[constraints[i]]=std::pair<int,T>(call_count,forces[i]);
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Increment_Forces(const Matrix<T,Dynamic,1>& forces)
{
    stored_forces.resize(forces.size());
    for(int i=0;i<constraints.size();i++){
        if(force_memory.count(constraints[i])){// this could be compressed if I could be sure that the force would be initialized properly
            auto& memory=force_memory[constraints[i]];
            memory.first=call_count;
            memory.second+=forces[i];
            stored_forces[i]=memory.second;
        }
        else{
            force_memory[constraints[i]]=std::pair<int,T>(call_count,forces[i]);
            stored_forces[i]=forces[i];
        }
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.Find("RIGID_STRUCTURE_DATA"));
    typedef Matrix<T,1,RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE> CONSTRAINT_VECTOR;
    std::vector<Triplet<CONSTRAINT_VECTOR>> terms;
    std::vector<T> rhs;
    RIGID_STRUCTURE_INDEX_MAP<TV> index_map;
    constraints.clear();
    for(int s1=0;s1<rigid_data->structures.size();s1++){
        for(int s2=s1+1;s2<rigid_data->structures.size();s2++){
            auto structure1=rigid_data->structures[s1];
            auto structure2=rigid_data->structures[s2];
            TV offset1,offset2;
            TV direction=structure1->Displacement(data,structure2,offset1,offset2);
            T distance=direction.norm();
            T constraint_violation=distance-structure1->collision_radius-structure2->collision_radius;
            T distance_condition=-.001;
            std::pair<int,T> remembered=force_memory[CONSTRAINT(s1,s2)];
            //std::cout<<s1<<" "<<rigid_structure1->name<<" "<<s2<<" "<<rigid_structure2->name<<": constraint violation "<<constraint_violation<<std::endl;
            //std::cout<<"s1 com: "<<rigid_structure1->frame.position.transpose()<<" s2 com: "<<rigid_structure2->frame.position.transpose()<<" r1:
            //"<<rigid_structure1->collision_radius<<" r2: "<<rigid_structure2->collision_radius<<std::endl;
            std::cout<<"Constraint violation: "<<constraint_violation<<" remembered call count: "<<remembered.first<<" call count: "<<call_count<<" remembered force: "<<remembered.second<<std::endl;
            if(constraint_violation<distance_condition || (remembered.first==call_count && remembered.second>0)){
                TV x1=structure1->frame.position+offset1;
                TV x2=structure2->frame.position+offset2;
                TV object_offset1=structure1->frame.orientation.inverse()*offset1;
                TV object_offset2=structure1->frame.orientation.inverse()*offset2;
                std::cout<<"Constraint between "<<s1<<" and "<<s2<<std::endl;
                T factor=500;
                CONSTRAINT_VECTOR DC_DA2=factor*RIGID_STRUCTURE_INDEX_MAP<TV>::DC_DA(*structure2,object_offset2,x1,x2,direction);
                CONSTRAINT_VECTOR DC_DA1=factor*RIGID_STRUCTURE_INDEX_MAP<TV>::DC_DA(*structure1,object_offset1,x1,x2,direction);
                terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s2,DC_DA2));
                terms.push_back(Triplet<CONSTRAINT_VECTOR>(constraints.size(),s1,-DC_DA1));
                Matrix<T,t+d,t+d> force_balance_contribution2=factor*remembered.second*RIGID_STRUCTURE_INDEX_MAP<TV>::DF_DA(*structure2,object_offset2,x1,x2,direction);
                Matrix<T,t+d,t+d> force_balance_contribution1=-factor*remembered.second*RIGID_STRUCTURE_INDEX_MAP<TV>::DF_DA(*structure1,object_offset1,x1,x2,direction);
                //Matrix<T,d,t+d> force_balance_contribution2=factor*remembered.second*RIGID_STRUCTURE_INDEX_MAP<TV>::DD_DV(*structure2,object_offset2,x1,x2,direction);
                //Matrix<T,d,t+d> force_balance_contribution1=-factor*remembered.second*RIGID_STRUCTURE_INDEX_MAP<TV>::DD_DV(*structure1,object_offset1,x1,x2,direction);
                std::cout<<"DRDA VOL: "<<force_balance_contribution1<<std::endl;
                  std::cout<<"END DRDA"<<std::endl;
                std::cout<<"F_i: "<<constraint_violation<<std::endl;
                auto hessian=factor*constraint_violation*RIGID_STRUCTURE_INDEX_MAP<TV>::Full_Hessian(structure1->twist.angular,structure2->twist.angular,offset1,offset2,x1,x2,direction);
                //std::cout<<"Hessian: "<<std::endl<<hessian<<std::endl;

                for(int j=0;j<t+d;j++){
                    for(int k=0;k<t+d;k++){
                        if(abs(force_balance_contribution1(j,k))>1e-6){
                            force_terms.push_back(Triplet<T>(s1*(t+d)+j,s1*(t+d)+k,force_balance_contribution1(j,k)));
                        }
                        if(abs(force_balance_contribution2(j,k))>1e-6){
                            force_terms.push_back(Triplet<T>(s2*(t+d)+j,s2*(t+d)+k,force_balance_contribution2(j,k)));
                        }
                        // Hessian terms
                        hessian_terms.push_back(Triplet<T>(s1*(t+d)+j,s1*(t+d)+k,hessian(j,k)));
                        hessian_terms.push_back(Triplet<T>(s1*(t+d)+j,s2*(t+d)+k,hessian(j,t+d+k)));
                        hessian_terms.push_back(Triplet<T>(s2*(t+d)+j,s1*(t+d)+k,hessian(t+d+j,k)));
                        hessian_terms.push_back(Triplet<T>(s2*(t+d)+j,s2*(t+d)+k,hessian(t+d+j,t+d+k)));
                    }
                }
                rhs.push_back(-factor*constraint_violation);
                std::cout<<"VOL RHS contribution: "<<(DC_DA1.transpose()*remembered.second).transpose()<<std::endl;
                right_hand_side.template block<t+d,1>(s1*(t+d),0)+=DC_DA1.transpose()*remembered.second;
                right_hand_side.template block<t+d,1>(s2*(t+d),0)-=DC_DA2.transpose()*remembered.second;
                constraints.push_back(CONSTRAINT(s1,s2));
            }
        }
    }
    stored_forces.resize(constraints.size(),1);
    for(int i=0;i<constraints.size();i++){
        stored_forces(i,0)=force_memory[constraints[i]].second;
    }
    constraint_terms.resize(constraints.size(),RIGID_STRUCTURE_INDEX_MAP<TV>::STATIC_SIZE*rigid_data->structures.size());
    constraint_rhs.resize(rhs.size(),1);
    for(int i=0;i<rhs.size();i++){constraint_rhs(i,0)=rhs[i];}
    Flatten_Matrix(terms,constraint_terms);
    call_count++;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void VOLUME_EXCLUSION_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
    osg::Group* group=node->asGroup();
    group->removeChild(getNamedChild(group,Static_Name()));
    auto rigid_data=std::static_pointer_cast<RIGID_STRUCTURE_DATA<TV>>(data.Find("RIGID_STRUCTURE_DATA"));
    osg::Group* volume_exclusion_group=new osg::Group();
    volume_exclusion_group->setName(Static_Name());
    std::cout<<"Volume exclusion constraints: "<<constraints.size()<<std::endl;
    for(int i=0;i<constraints.size();i++){
        auto lineGeometry=new osg::Geometry();
        auto vertices=new osg::Vec3Array(2);
        const CONSTRAINT& constraint=constraints[i];
        int body_index1=constraint.first;
        int body_index2=constraint.second;
        std::cout<<"Between "<<body_index1<<" and "<<body_index2<<std::endl;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        auto firstAttachment=rigid_structure1->frame.position;
        auto secondAttachment=rigid_structure2->frame.position;
        (*vertices)[0].set(firstAttachment(0),firstAttachment(1),firstAttachment(2));
        (*vertices)[1].set(secondAttachment(0),secondAttachment(1),secondAttachment(2));
        lineGeometry->setVertexArray(vertices);
        auto colors=new osg::Vec4Array;
        colors->push_back(osg::Vec4(0.0f,0.0f,1.0f,1.0f));
        lineGeometry->setColorArray(colors);
        lineGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
        auto normals=new osg::Vec3Array;
        normals->push_back(osg::Vec3f(0.0f,-1.0f,0.0f));
        lineGeometry->setNormalArray(normals);
        lineGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);
        
        lineGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,2));
        auto lineGeode=new osg::Geode();
        lineGeode->addDrawable(lineGeometry);
        volume_exclusion_group->addChild(lineGeode);
    }
    group->addChild(volume_exclusion_group);
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(VOLUME_EXCLUSION_CONSTRAINT)
GENERIC_TYPE_DEFINITION(VOLUME_EXCLUSION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(VOLUME_EXCLUSION_CONSTRAINT,void)
{
    auto volume_exclusion_constraint=std::make_shared<VOLUME_EXCLUSION_CONSTRAINT<TV>>();
    simulation.force.push_back(volume_exclusion_constraint);
    return 0;
}
