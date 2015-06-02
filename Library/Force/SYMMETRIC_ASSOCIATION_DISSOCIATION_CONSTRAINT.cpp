#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Force/SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
#include <Force/FORCE.h>
#include <Indexing/RIGID_STRUCTURE_INDEX_MAP.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/OSG_HELPERS.h>
#include <Utilities/RANDOM.h>
#include <unsupported/Eigen/BVH>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/LineWidth>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> std::shared_ptr<FORCE_REFERENCE<typename TV::Scalar>> SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Create_Stored_Force() const
{
    return std::static_pointer_cast<FORCE_REFERENCE<T>>(std::make_shared<STORED_SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<STORED_SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    information->constraints=constraints;
    information->value.resize(information->Size());
    information->value.setZero();
    for(int i=0;i<information->constraints.size();i++){
        FORCE_VECTOR& value=force_memory[information->constraints[i]].second;
        information->value.template block<d,1>(d*i,0)=value.template block<d,1>(0,0);;
        information->value.template block<t,1>(d*information->constraints.size()+t*i,0)=value.template block<t,1>(d,0);;
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<const STORED_SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        force_memory[information->constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,information->value.template block<d+t,1>(i*(d+t),0));
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment)
{
    call_count+=increment;
    auto information=std::static_pointer_cast<STORED_SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        FORCE_VECTOR value;
        value.template block<d,1>(0,0)=information->value.template block<d,1>(d*i,0);
        value.template block<t,1>(d,0)=information->value.template block<t,1>(d*information->constraints.size()+t*i,0);
        if(force_memory.find(information->constraints[i])!=force_memory.end()){
            auto& memory=force_memory[information->constraints[i]];
            memory.first=call_count;
            memory.second+=increment*value;
        }
        else{
            force_memory[information->constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,increment*value);
        }
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> ROTATION<TV> SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Find_Appropriate_Rotation(const ROTATION<TV>& rotation1,const ROTATION<TV>& rotation2)
{
    return ROTATION<TV>(rotation1.inverse()*(rotation2*rotation1.inverse()).inverse()).Scale_Angle((T).5);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    if(stochastic){
        for(auto memory : force_memory){memory.second.second.setZero();}
        constraints.clear();

        for(int i=0;i<interaction_types.size();i++){
            // build acceleration structure out of all binding sites
            auto& interaction_type=interaction_types[i];
            auto& sites=interaction_type.sites;
            std::vector<AlignedBox<T,d>> bounding_list(sites.size());
            std::vector<int> index_list(sites.size());
            for(int s=0;s<sites.size();s++){
                auto structure=rigid_data->structures[std::get<0>(sites[s])];
                bounding_list[s]=AlignedBox<T,d>(structure->frame*std::get<1>(sites[s])-TV::Constant(interaction_type.bond_distance_threshold),
                    structure->frame*std::get<1>(sites[s])+TV::Constant(interaction_type.bond_distance_threshold));
                index_list[s]=s;
            }
            KdBVH<T,3,int> hierarchy_second_site(index_list.begin(),index_list.end(),bounding_list.begin(),bounding_list.end());

            for(int candidate_first=0;candidate_first<interaction_type.sites.size();candidate_first++){
                auto first_site=interaction_type.sites[candidate_first];
                int s1=std::get<0>(first_site);
                auto structure1=rigid_data->structures[s1];
                auto first_site_position=structure1->frame*std::get<1>(first_site);
                PROXIMITY_SEARCH<TV> proximity_search(data,first_site_position,interaction_type.bond_distance_threshold);
                ROTATION<TV> binder1_frame=structure1->frame.orientation;
                BVIntersect(hierarchy_second_site,proximity_search);
                LOG::cout<<proximity_search.candidates.size()<<" candidates"<<std::endl;
                for(auto candidate_second : proximity_search.candidates){
                    auto second_site=interaction_type.sites[candidate_second];
                    int s2=std::get<0>(second_site);
                    if(s1==s2){continue;}
                    std::pair<int,int> partnership(std::min(s1,s2),std::max(s1,s2));
                    bool constraint_active=false;
                    CONSTRAINT constraint(i,candidate_first,candidate_second);
                    std::pair<int,FORCE_VECTOR> remembered=force_memory[constraint];
                    if(remembered.first==call_count){
                        T dissociation_rate=1/interaction_type.base_dissociation_time;
                        T cumulative_distribution=1-exp(-dissociation_rate*dt);
                        constraint_active=data.random.Uniform((T)0,(T)1)>cumulative_distribution;
                        if(!constraint_active){
                            std::get<2>(interaction_type.sites[candidate_first])=false;
                            std::get<2>(interaction_type.sites[candidate_second])=false;
                            std::cout<<"CONSTRAINT DEACTIVATED: "<<s1<<" "<<s2<<std::endl;
                            partners[partnership]=false;
                            partner_interactions.erase(partnership);
                        }
                    }
                    else if(!std::get<2>(interaction_type.sites[candidate_first]) && !std::get<2>(second_site) && !partners[partnership]){ // do not re-bind already bound
                        auto structure2=rigid_data->structures[s2];
                        auto second_site_position=structure2->frame*std::get<1>(second_site);
                        ROTATION<TV> binder2_frame=structure2->frame.orientation;
                        LOG::cout<<"First axis: "<<structure2->frame.orientation.Axis().transpose()<<" angle: "<<structure2->frame.orientation.Angle()<<std::endl;
                        LOG::cout<<"First site: "<<first_site_position.transpose()<<" Second site: "<<second_site_position.transpose()<<std::endl;
                        TV direction=data.Minimum_Offset(first_site_position,second_site_position);
                        LOG::cout<<"Direction: "<<direction.transpose()<<std::endl;
                        T bond_distance=direction.norm();
                        T orientation_compatibility=T(),position_compatibility=T();
                        if(bond_distance<interaction_type.bond_distance_threshold){
                            position_compatibility=1-bond_distance/interaction_type.bond_distance_threshold;
                            ROTATION<TV> composed_rotation(binder2_frame.inverse()*binder1_frame*interaction_type.relative_orientation.inverse());
                            ROTATION<TV> relative_rotation(binder2_frame.inverse()*binder1_frame);
                            LOG::cout<<"Angles: "<<structure2->frame.orientation.Angle()<<" "<<structure1->frame.orientation.Angle()<<" "<<interaction_type.relative_orientation.Angle()<<" "<<(structure2->frame.orientation.inverse()*structure1->frame.orientation*interaction_type.relative_orientation.inverse()).w()<<" "<<(structure2->frame.orientation.inverse()*structure1->frame.orientation).w()<<std::endl;
                            LOG::cout<<"Relative rotation axis: "<<relative_rotation.Axis().transpose()<<" Angle: "<<relative_rotation.Angle()<<std::endl;;
                            LOG::cout<<"Desired relative rotation axis: "<<interaction_type.relative_orientation.Axis().transpose()<<" Angle: "<<interaction_type.relative_orientation.Angle()<<std::endl;
                            LOG::cout<<"Composed angle: "<<composed_rotation.Angle()<<std::endl;
                            ROTATION<TV> R2_O=binder2_frame*interaction_type.relative_orientation;
                            LOG::cout<<"R2_0: Axis: "<<R2_O.Axis().transpose()<<" angle: "<<R2_O.Angle()<<std::endl;
                            LOG::cout<<"R1: Axis: "<<binder1_frame.Axis().transpose()<<" angle: "<<binder1_frame.Angle()<<std::endl;
                            T angle_magnitude=std::abs(composed_rotation.Angle());
                            orientation_compatibility=std::max((T)0,1-std::min(angle_magnitude,2*(T)M_PI-angle_magnitude)/interaction_type.bond_orientation_threshold);}
                        T compatibility=orientation_compatibility*position_compatibility;
                        T association_rate=compatibility/interaction_type.base_association_time;
                        T cumulative_distribution=1-exp(-association_rate*dt);
                        constraint_active=data.random.Uniform((T)0,(T)1)<cumulative_distribution;
                        LOG::cout<<"Maybe activating constraint: "<<constraint_active<<" compatibility "<<compatibility<<" bond_distance: "<<bond_distance<<" orientation_compatibility: "<<orientation_compatibility<<std::endl;
                        if(constraint_active){
                            std::cout<<"CONSTRAINT ACTIVATED: "<<s1<<" "<<s2<<" "<<candidate_first<<" "<<candidate_second<<" remembered: "<<remembered.first<<" call count: "<<call_count<<" interaction: "<<i<<" test: "<<std::get<2>(interaction_type.sites[candidate_first])<<" "<<std::get<2>(interaction_type.sites[candidate_second])<<" "<<partners[partnership]<<std::endl;
                            std::get<2>(interaction_type.sites[candidate_first])=true;
                            std::get<2>(interaction_type.sites[candidate_second])=true;
                            std::cout<<"CONSTRAINT ACTIVATED: "<<s1<<" "<<s2<<" "<<candidate_first<<" "<<candidate_second<<" remembered: "<<remembered.first<<" call count: "<<call_count<<" interaction: "<<i<<" test: "<<std::get<2>(interaction_type.sites[candidate_first])<<" "<<std::get<2>(interaction_type.sites[candidate_second])<<std::endl;
                            partners[partnership]=true;
                            partner_interactions[partnership]=i;
                            std::cout<<"ACTIVE interactions:"<<std::endl;
                            for(auto active_pair : partner_interactions){
                                std::cout<<"pair: "<<active_pair.first.first<<", "<<active_pair.first.second<<" interaction: "<<active_pair.second<<std::endl;
                            }
                        }
                    }
                    if(constraint_active){constraints.push_back(constraint);}}
            }
        }}

    Matrix<T,t,t+d> angular_to_constraint;angular_to_constraint.setZero();
    angular_to_constraint.template block<t,t>(0,d).setIdentity();
    std::vector<Triplet<LINEAR_CONSTRAINT_MATRIX>> linear_terms;
    std::vector<Triplet<ANGULAR_CONSTRAINT_MATRIX>> angular_terms;
    constraint_rhs.resize(constraints.size()*(d+t));
    stored_forces.resize(constraints.size()*(d+t));
    for(int i=0;i<constraints.size();i++){
        auto interaction_index=constraints[i];
        auto interaction_type=interaction_types[std::get<0>(interaction_index)];
        int first_site_index=std::get<1>(interaction_index);
        int second_site_index=std::get<2>(interaction_index);
        auto first_site=interaction_type.sites[first_site_index];
        auto second_site=interaction_type.sites[second_site_index];
        auto s1=std::get<0>(first_site);
        auto s2=std::get<0>(second_site);
        auto v1=std::get<1>(first_site);
        auto v2=std::get<1>(second_site);
        auto structure1=rigid_data->structures[s1];
        auto structure2=rigid_data->structures[s2];
        FORCE_VECTOR rhs;
        //std::cout<<"Site indices: "<<s1<<" "<<s2<<" "<<first_site_index<<" "<<second_site_index<<std::endl;
        
        auto x1=structure1->frame*v1;
        auto x2=structure2->frame*v2;
        //TV direction=structure1->Displacement(data,*structure2,offset1,offset2).normalized(); // use core-core direction for stability reasons
        TV direction=data.Minimum_Offset(x1,x2); // can't easily use point-point distance then, though.
        
        LINEAR_CONSTRAINT_MATRIX dC_dX1=RIGID_STRUCTURE_INDEX_MAP<TV>::Velocity_Map(*structure1,v1);
        LINEAR_CONSTRAINT_MATRIX dC_dX2=RIGID_STRUCTURE_INDEX_MAP<TV>::Velocity_Map(*structure2,v2);
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,s2,dC_dX2));
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,s1,-dC_dX1));
        rhs.template block<d,1>(0,0)=-direction;

        ROTATION<TV> relative_orientation=interaction_type.relative_orientation;

        ROTATION<TV> R1_current=ROTATION<TV>::From_Rotation_Vector(structure1->twist.angular);
        ROTATION<TV> R2_current=ROTATION<TV>::From_Rotation_Vector(structure2->twist.angular);
        ROTATION<TV> R1_base=(R1_current.inverse()*structure1->frame.orientation)*relative_orientation.inverse();
        ROTATION<TV> R2_base=R2_current.inverse()*structure2->frame.orientation;
        ROTATION<TV> RC=Find_Appropriate_Rotation(R1_current*R1_base,R2_current*R2_base);
        T_SPIN first_rotation_error_vector,second_rotation_error_vector;
        
        // R1=R1_s*R1_0
        // R1_s*O*R1_0

        // R1
        ANGULAR_CONSTRAINT_MATRIX dC_dA1=RIGID_STRUCTURE_INDEX_MAP<TV>::Construct_Constraint_Matrix(R1_current,R1_base*RC,first_rotation_error_vector)*angular_to_constraint;
        // R2*O
        ANGULAR_CONSTRAINT_MATRIX dC_dA2=RIGID_STRUCTURE_INDEX_MAP<TV>::Construct_Constraint_Matrix(R2_current,R2_base*RC,second_rotation_error_vector)*angular_to_constraint;
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,s2,dC_dA2.eval()));
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,s1,-dC_dA1.eval()));
        //R2^-1*R1*O^-1=I
        //R1*O^-1=R2
        //R1=R2*O
        ROTATION<TV> composed_rotation(relative_orientation.inverse()*(structure2->frame.orientation*RC).inverse()*(structure1->frame.orientation*RC));
        LOG::cout<<"Composed angle: "<<composed_rotation.Angle()<<" axis: "<<composed_rotation.Axis().transpose()<<std::endl;
        LOG::cout<<"Composed angle w: "<<composed_rotation.Vec().transpose()<<" W: "<<composed_rotation.W()<<std::endl;
        T_SPIN total_rotation_error=second_rotation_error_vector-first_rotation_error_vector;
        LOG::cout<<"rotation error: "<<total_rotation_error.transpose()<<std::endl;
        rhs.template block<t,1>(d,0)=-total_rotation_error;
        constraint_rhs.template block<d,1>(d*i,0)=rhs.template block<d,1>(0,0);
        constraint_rhs.template block<t,1>(d*constraints.size()+i*t,0)=rhs.template block<t,1>(d,0);
        //constraint_rhs.template block<d+t,1>((d+t)*i,0)=rhs;
        auto& remembered=force_memory[interaction_index];
        if(remembered.first!=call_count){remembered.second.setZero();}
        stored_forces.template block<d+t,1>((d+t)*i,0)=remembered.second;
        TV right_hand_force=remembered.second.template block<d,1>(0,0);
        T_SPIN right_hand_torque=remembered.second.template block<t,1>(d,0);
        right_hand_side.template block<d+t,1>(s1*(d+t),0)+=dC_dA1.transpose()*right_hand_torque+dC_dX1.transpose()*right_hand_force;
        right_hand_side.template block<d+t,1>(s2*(d+t),0)-=dC_dA2.transpose()*right_hand_torque+dC_dX2.transpose()*right_hand_force;
    }
    constraint_terms.resize(constraints.size()*(d+t),rigid_data->Velocity_DOF());
    //Flattern_Matrix(linear_terms,constraint_terms);
    Flatten_Matrices(linear_terms,d*constraints.size(),angular_terms,constraint_terms);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
    std::cout<<"Viewer has "<<interaction_types.size()<<" types "<<constraints.size()<<" constraints"<<std::endl;
    osg::Group* group=node->asGroup();
    group->removeChild(getNamedChild(group,Static_Name()));
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    osg::Group* volume_exclusion_group=new osg::Group();
    volume_exclusion_group->setName(Static_Name());
    for(int i=0;i<constraints.size();i++){
        auto lineGeometry=new osg::Geometry();
        auto vertices=new osg::Vec3Array(2);

        auto interaction_index=constraints[i];
        auto interaction_type=interaction_types[std::get<0>(interaction_index)];
        int first_site_index=std::get<1>(interaction_index);
        int second_site_index=std::get<2>(interaction_index);
        auto first_site=interaction_type.sites[first_site_index];
        auto second_site=interaction_type.sites[second_site_index];
        auto body_index1=std::get<0>(first_site);
        auto body_index2=std::get<0>(second_site);
        auto v1=std::get<1>(first_site);
        auto v2=std::get<1>(second_site);

        std::cout<<"Between "<<body_index1<<" and "<<body_index2<<std::endl;
        auto rigid_structure1=rigid_data->structures[body_index1];
        auto rigid_structure2=rigid_data->structures[body_index2];
        auto firstAttachment=rigid_structure1->frame*(v1/2);
        auto secondAttachment=firstAttachment+data.Minimum_Offset(firstAttachment,rigid_structure2->frame*(v2/2));
        //auto firstAttachment=rigid_structure1->frame.position;
        //auto secondAttachment=firstAttachment+data.Minimum_Offset(firstAttachment,rigid_structure2->frame.position);
        (*vertices)[0].set(firstAttachment(0),firstAttachment(1),firstAttachment(2));
        (*vertices)[1].set(secondAttachment(0),secondAttachment(1),secondAttachment(2));
        lineGeometry->setVertexArray(vertices);
        auto colors=new osg::Vec4Array;
        osg::Vec4 color;
        color[3]=1;
        color[std::get<0>(interaction_index)]=1;
        colors->push_back(color);
        //colors->push_back(osg::Vec4(1.0f,0.0f,0.0f,1.0f));
        lineGeometry->setColorArray(colors);
        lineGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
        auto normals=new osg::Vec3Array;
        normals->push_back(osg::Vec3f(0.0f,-1.0f,0.0f));
        lineGeometry->setNormalArray(normals);
        lineGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);
        osg::StateSet* stateset=new osg::StateSet;
        osg::LineWidth* lineWidth=new osg::LineWidth();
        lineWidth->setWidth(4.0f);
        stateset->setAttributeAndModes(lineWidth,osg::StateAttribute::ON);
        stateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
        lineGeometry->setStateSet(stateset);
        lineGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,2));
        auto lineGeode=new osg::Geode();
        lineGeode->addDrawable(lineGeometry);
        volume_exclusion_group->addChild(lineGeode);
    }
    group->addChild(volume_exclusion_group);
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT,void)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto constraint=simulation.force.template Find_Or_Create<SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>();
    auto interactions=node["interactions"];
    for(Json::ValueIterator it=interactions.begin();it!=interactions.end();it++){
        typename SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::INTERACTION_TYPE interaction;
        interaction.bond_distance_threshold=(*it)["bond_distance_threshold"].asDouble();
        interaction.bond_orientation_threshold=(*it)["bond_orientation_threshold"].asDouble();
        interaction.base_association_time=(*it)["base_association_time"].asDouble();
        interaction.base_dissociation_time=(*it)["base_dissociation_time"].asDouble();
        TV site_offset;Parse_Vector((*it)["site_offset"],site_offset);
        {
            auto& binder_orientation=(*it)["binder_orientation"];
            TV primary_axis;Parse_Vector(binder_orientation["primary_axis"],primary_axis);
            TV secondary_axis_from;Parse_Vector(binder_orientation["secondary_axis"]["from"],secondary_axis_from);
            TV secondary_axis_to;Parse_Vector(binder_orientation["secondary_axis"]["to"],secondary_axis_to);
            ROTATION<TV> primary_rotation=ROTATION<TV>::From_Rotated_Vector(site_offset,primary_axis);
            ROTATION<TV> secondary_rotation=ROTATION<TV>::From_Rotated_Vector_Around_Axis(primary_rotation*secondary_axis_from,secondary_axis_to,primary_axis);
            interaction.relative_orientation=secondary_rotation*primary_rotation;
        }
        auto sites=(*it)["sites"];
        for(auto site_it=sites.begin();site_it!=sites.end();site_it++){
            std::tuple<int,TV,bool> site;
            std::get<0>(site)=rigid_data->Structure_Index((*site_it)["name"].asString());
            std::get<1>(site)=site_offset;
            std::get<2>(site)=false;
            interaction.sites.push_back(site);
        }
        constraint->interaction_types.push_back(interaction);
    }
    return 0;
}
