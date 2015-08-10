#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Force/ASSOCIATION_DISSOCIATION_CONSTRAINT.h>
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
template<class TV> std::shared_ptr<FORCE_REFERENCE<typename TV::Scalar>> ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Create_Stored_Force() const
{
    return std::static_pointer_cast<FORCE_REFERENCE<T>>(std::make_shared<STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Pack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    information->constraints=constraints;
    information->value.resize(information->Size());
    information->value.setZero();
    for(int i=0;i<information->constraints.size();i++){
        FORCE_VECTOR& value=force_memory[information->constraints[i]].second;
        information->value.template block<d,1>(d*i,0)=value.template block<d,1>(0,0);;
        information->value.template block<t,1>(d*information->constraints.size()+t*i,0)=value.template block<t,1>(d,0);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Unpack_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information)
{
    auto information=std::static_pointer_cast<const STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        force_memory[information->constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,information->value.template block<d+t,1>(i*(d+t),0));}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Increment_Forces(std::shared_ptr<FORCE_REFERENCE<T>> force_information,int increment)
{
    call_count+=increment;
    auto information=std::static_pointer_cast<STORED_ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>(force_information);
    for(int i=0;i<information->constraints.size();i++){
        FORCE_VECTOR value;
        value.template block<d,1>(0,0)=information->value.template block<d,1>(d*i,0);
        value.template block<t,1>(d,0)=information->value.template block<t,1>(d*information->constraints.size()+t*i,0);
        if(force_memory.find(information->constraints[i])!=force_memory.end()){
            auto& memory=force_memory[information->constraints[i]];
            memory.first=call_count;
            memory.second+=increment*value;}
        else{force_memory[information->constraints[i]]=std::pair<int,FORCE_VECTOR>(call_count,increment*value);}}
}
///////////////////////////////////////////////////////////////////////
template<class TV> ROTATION<TV> ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Find_Appropriate_Rotation(const ROTATION<TV>& rotation1,const ROTATION<TV>& rotation2)
{
    return ROTATION<TV>(rotation1.inverse()*(rotation2*rotation1.inverse()).inverse()).Scale_Angle((T).5);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Interaction_Candidates(DATA<TV>& data,const T dt,int type_index,int type1,int type2)
{
    auto& interaction_type=interaction_types[type_index];
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    // build acceleration structure out of all binding sites
    auto& first_sites=interaction_type.sites[type1];
    auto& second_sites=interaction_type.sites[type2];
    std::vector<AlignedBox<T,d>> bounding_list(second_sites.size());
    std::vector<int> index_list(second_sites.size());
    for(int s=0;s<second_sites.size();s++){
        TV site_location=rigid_data->structures[std::get<0>(second_sites[s])]->frame*interaction_type.site_offsets[type2];
        bounding_list[s]=AlignedBox<T,d>(site_location,site_location);
        index_list[s]=s;}
    KdBVH<T,3,int> hierarchy_second_site(index_list.begin(),index_list.end(),bounding_list.begin(),bounding_list.end());

    for(int candidate_first=0;candidate_first<first_sites.size();candidate_first++){
        auto& first_site=first_sites[candidate_first];
        int s1=std::get<0>(first_site);
        auto structure1_frame=rigid_data->structures[s1]->frame;
        auto first_site_position=structure1_frame*interaction_type.site_offsets[type1];
        PROXIMITY_SEARCH<TV> proximity_search(data,first_site_position,interaction_type.bond_distance_threshold);
        BVIntersect(hierarchy_second_site,proximity_search);
        for(auto candidate_second : proximity_search.candidates){
            auto& second_site=second_sites[candidate_second];
            int s2=std::get<0>(second_site);
            // prevent a) two bindings between a pair of bodies b) double-counting in the symmetric site case
            if(s1==s2 || (type1==type2 && s1>s2)){continue;}
            std::pair<int,int> partnership(std::min(s1,s2),std::max(s1,s2));
            bool constraint_active=false;
            CONSTRAINT constraint(type_index,candidate_first,candidate_second);
            std::pair<int,FORCE_VECTOR> remembered=force_memory[constraint];
            if(remembered.first==call_count){
                T dissociation_rate=1/interaction_type.base_dissociation_time;
                T cumulative_distribution=1-exp(-dissociation_rate*dt);
                constraint_active=data.random.Uniform((T)0,(T)1)>cumulative_distribution;
                if(!constraint_active){
                    std::get<1>(first_site)=false;
                    std::get<1>(second_site)=false;
                    partners[partnership]=false;
                    std::cout<<"CONSTRAINT DEACTIVATED: "<<s1<<" "<<s2<<std::endl;}}
            else if(!std::get<1>(first_site) && !std::get<1>(second_site) && !partners[partnership]){ // do not re-bind already bound
                auto structure2_frame=rigid_data->structures[s2]->frame;
                auto second_site_position=structure2_frame*interaction_type.site_offsets[type2];
                T bond_distance=data.Minimum_Offset(first_site_position,second_site_position).norm();
                T orientation_compatibility=T(),position_compatibility=T();
                if(bond_distance<interaction_type.bond_distance_threshold){
                    position_compatibility=1-bond_distance/interaction_type.bond_distance_threshold;
                    ROTATION<TV> composed_rotation(structure2_frame.orientation.inverse()*structure1_frame.orientation*interaction_type.relative_orientation.inverse());
                    T angle_magnitude=std::abs(composed_rotation.Angle());
                    orientation_compatibility=std::max((T)0,1-std::min(angle_magnitude,2*(T)M_PI-angle_magnitude)/interaction_type.bond_orientation_threshold);
                }
                T association_rate=orientation_compatibility*position_compatibility/interaction_type.base_association_time;
                T cumulative_distribution=1-exp(-association_rate*dt);
                constraint_active=data.random.Uniform((T)0,(T)1)<cumulative_distribution;
                if(constraint_active){
                    std::cout<<"CONSTRAINT ACTIVATED: "<<s1<<" "<<s2<<std::endl;
                    //std::cout<<"CONSTRAINT ACTIVATED: "<<s1<<" "<<s2<<" "<<candidate_first<<" "<<candidate_second<<" remembered: "<<remembered.first<<" call count: "<<call_count<<" interaction: "<<i<<" test: "<<std::get<1>(interaction_type.sites[candidate_first])<<" "<<std::get<1>(interaction_type.sites[candidate_second])<<" "<<partners[partnership]<<std::endl;
                    std::get<1>(first_site)=true;
                    std::get<1>(second_site)=true;
                    partners[partnership]=true;}}
            if(constraint_active){constraints.push_back(constraint);}}}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Linearize(DATA<TV>& data,const T dt,const T target_time,std::vector<Triplet<T>>& hessian_terms,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& constraint_terms,SparseMatrix<T>& constraint_forces,Matrix<T,Dynamic,1>& right_hand_side,Matrix<T,Dynamic,1>& constraint_rhs,bool stochastic)
{
    auto rigid_data=data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    if(stochastic){
        for(auto memory : force_memory){memory.second.second.setZero();}
        constraints.clear();
        for(int i=0;i<interaction_types.size();i++){
            int site_count=interaction_types[i].sites.size();
            if(site_count==1){Interaction_Candidates(data,dt,i,0,0);}
            else{Interaction_Candidates(data,dt,i,0,1);}}}

    Matrix<T,t,t+d> angular_to_constraint;angular_to_constraint.setZero();
    angular_to_constraint.template block<t,t>(0,d).setIdentity();
    std::vector<Triplet<LINEAR_CONSTRAINT_MATRIX>> linear_terms;
    std::vector<Triplet<ANGULAR_CONSTRAINT_MATRIX>> angular_terms;
    constraint_rhs.resize(constraints.size()*(d+t));
    stored_forces.resize(constraints.size()*(d+t));
    for(int i=0;i<constraints.size();i++){
        auto interaction_index=constraints[i];
        auto interaction_type=interaction_types[std::get<0>(interaction_index)];
        int first_site_index=std::get<1>(interaction_index),second_site_index=std::get<2>(interaction_index);
        auto first_site=interaction_type.sites[0][first_site_index],second_site=interaction_type.sites.back()[second_site_index];
        auto s1=std::get<0>(first_site),s2=std::get<0>(second_site);
        auto v1=interaction_type.site_offsets[0],v2=interaction_type.site_offsets.back();
        auto structure1=rigid_data->structures[s1],structure2=rigid_data->structures[s2];
        
        //TV direction=structure1->Displacement(data,*structure2,offset1,offset2).normalized(); // use core-core direction for stability reasons
        TV direction=data.Minimum_Offset(structure1->frame*v1,structure2->frame*v2); // can't easily use point-point distance then, though.
        
        LINEAR_CONSTRAINT_MATRIX dC_dX1=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(*structure1,v1);
        LINEAR_CONSTRAINT_MATRIX dC_dX2=RIGID_STRUCTURE_INDEX_MAP<TV>::Map_Twist_To_Velocity(*structure2,v2);
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,s2,dC_dX2));
        linear_terms.push_back(Triplet<LINEAR_CONSTRAINT_MATRIX>(i,s1,-dC_dX1));

        ROTATION<TV> relative_orientation_inverse=interaction_type.relative_orientation.inverse();
        ROTATION<TV> R1_current=ROTATION<TV>::From_Rotation_Vector(structure1->twist.angular);
        ROTATION<TV> R2_current=ROTATION<TV>::From_Rotation_Vector(structure2->twist.angular);
        ROTATION<TV> R1_base=(R1_current.inverse()*structure1->frame.orientation)*relative_orientation_inverse;
        ROTATION<TV> R2_base=R2_current.inverse()*structure2->frame.orientation;
        ROTATION<TV> RC=Find_Appropriate_Rotation(R1_current*R1_base,R2_current*R2_base);
        T_SPIN first_rotation_error_vector,second_rotation_error_vector;
        
        ANGULAR_CONSTRAINT_MATRIX dC_dA1=RIGID_STRUCTURE_INDEX_MAP<TV>::Orientation_Constraint_Matrix(R1_current,R1_base*RC,first_rotation_error_vector)*angular_to_constraint;
        ANGULAR_CONSTRAINT_MATRIX dC_dA2=RIGID_STRUCTURE_INDEX_MAP<TV>::Orientation_Constraint_Matrix(R2_current,R2_base*RC,second_rotation_error_vector)*angular_to_constraint;
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,s2,dC_dA2.eval()));
        angular_terms.push_back(Triplet<ANGULAR_CONSTRAINT_MATRIX>(i,s1,-dC_dA1.eval()));
        ROTATION<TV> composed_rotation(relative_orientation_inverse*(structure2->frame.orientation*RC).inverse()*(structure1->frame.orientation*RC));
        T_SPIN total_rotation_error=second_rotation_error_vector-first_rotation_error_vector;
        constraint_rhs.template block<d,1>(d*i,0)=-direction;
        constraint_rhs.template block<t,1>(d*constraints.size()+i*t,0)=-total_rotation_error;
        auto& remembered=force_memory[interaction_index];
        if(remembered.first!=call_count){remembered.second.setZero();}
        stored_forces.template block<d+t,1>((d+t)*i,0)=remembered.second;
        TV right_hand_force=remembered.second.template block<d,1>(0,0);
        T_SPIN right_hand_torque=remembered.second.template block<t,1>(d,0);
        right_hand_side.template block<d+t,1>(s1*(d+t),0)+=dC_dA1.transpose()*right_hand_torque+dC_dX1.transpose()*right_hand_force;
        right_hand_side.template block<d+t,1>(s2*(d+t),0)-=dC_dA2.transpose()*right_hand_torque+dC_dX2.transpose()*right_hand_force;}
    constraint_terms.resize(constraints.size()*(d+t),rigid_data->Velocity_DOF());
    Flatten_Matrices(linear_terms,d*constraints.size(),angular_terms,constraint_terms);
}
///////////////////////////////////////////////////////////////////////
template<class TV> void ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::
Viewer(const DATA<TV>& data,osg::Node* node)
{
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
        int first_site_index=std::get<1>(interaction_index),second_site_index=std::get<2>(interaction_index);
        auto first_site=interaction_type.sites[0][first_site_index],second_site=interaction_type.sites.back()[second_site_index];
        int body_index1=std::get<0>(first_site),body_index2=std::get<0>(second_site);
        TV v1=interaction_type.site_offsets[0],v2=interaction_type.site_offsets.back();
        auto rigid_structure1=rigid_data->structures[body_index1],rigid_structure2=rigid_data->structures[body_index2];
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
GENERIC_TYPE_DEFINITION(ASSOCIATION_DISSOCIATION_CONSTRAINT)
DEFINE_AND_REGISTER_PARSER(ASSOCIATION_DISSOCIATION_CONSTRAINT,void)
{
    auto rigid_data=simulation.data.template Find<RIGID_STRUCTURE_DATA<TV>>();
    auto constraint=simulation.force.template Find_Or_Create<ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>>();
    auto interactions=node["interactions"];
    for(Json::ValueIterator it=interactions.begin();it!=interactions.end();it++){
        typename ASSOCIATION_DISSOCIATION_CONSTRAINT<TV>::INTERACTION_TYPE interaction;
        interaction.bond_distance_threshold=(*it)["bond_distance_threshold"].asDouble();
        interaction.bond_orientation_threshold=(*it)["bond_orientation_threshold"].asDouble();
        interaction.base_association_time=(*it)["base_association_time"].asDouble();
        interaction.base_dissociation_time=(*it)["base_dissociation_time"].asDouble();
        auto sites=(*it)["sites"];
        int site_types=0;
        for(auto site_it=sites.begin();site_it!=sites.end();site_it++,site_types++){
            interaction.sites.resize(site_types+1);
            interaction.site_offsets.resize(site_types+1);
            Parse_Vector((*site_it)["site_offset"],interaction.site_offsets[site_types]);
            auto bodies=(*site_it)["bodies"];
            for(auto body_it=bodies.begin();body_it!=bodies.end();body_it++){
                interaction.sites[site_types].push_back(std::pair<int,bool>(rigid_data->Structure_Index((*body_it).asString()),false));}}
        auto& binder_orientation=(*it)["binder_orientation"];
        TV primary_axis;Parse_Vector(binder_orientation["primary_axis"],primary_axis);
        TV secondary_axis_from;Parse_Vector(binder_orientation["secondary_axis"]["from"],secondary_axis_from);
        TV secondary_axis_to;Parse_Vector(binder_orientation["secondary_axis"]["to"],secondary_axis_to);
        ROTATION<TV> primary_rotation=ROTATION<TV>::From_Rotated_Vector(interaction.site_offsets[0],primary_axis);
        ROTATION<TV> secondary_rotation=ROTATION<TV>::From_Rotated_Vector_Around_Axis(primary_rotation*secondary_axis_from,secondary_axis_to,primary_axis);
        interaction.relative_orientation=secondary_rotation*primary_rotation;
        constraint->interaction_types.push_back(interaction);}
    return 0;
}
