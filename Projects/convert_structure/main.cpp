#include <Driver/DRIVER.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSE_SCENE.h>
#include <Utilities/LOG.h>
#include <algorithm>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <string>
#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/molecular_system.h>
#include "Vec3.h"
#include "OBB.h"
#include <typeinfo>
#include <string>
#include <algorithm>
//#include <json/json.h>
//#include <json/json-forwards.h>

/*#include <osg/ref_ptr>
#include <osgDB/Registry>
#include <osgDB/WriteFile>
#include <osg/Notify>*/
#include <iostream>
/*#include <osg/Geode>
#include <osg/Geometry>
#include <osgUtil/DelaunayTriangulator>*/

using namespace Mechanics;

// Defines a generic classifier that associates properties to objects. In this case, atom radius to atoms
typedef ESBTL::Generic_classifier<ESBTL::Radius_of_atom<double,ESBTL::Default_system::Atom> > T_Atom_classifier;

// defines an occupancy policy: discard atoms with no alternate location identification and with occupancy != 1
typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

// determines direction of bounding box rotation by applying forward and inverse rotation to all points and observing when points fall outside of aligned box
void verifyRotation(OBB obb_hull, std::vector<Vec3> points);

int main(int argc,char **argv)
{
    typedef double T;

    // defines a class (OBB_System) representing a molecular system
    typedef ESBTL::Molecular_system<ESBTL::Default_system_items,Vec3> OBB_System;
    
    // defines a class used to define what the system is made of (in this case two systems;
    // heavy atoms of water/nonwater molecules)
    typedef ESBTL::PDB_line_selector_two_systems Line_Selector;
    Line_Selector sel;
  
    std::vector<OBB_System> systems;

    //declare the build that will contruct the system from the pdb file.
    typedef ESBTL::All_atom_system_builder<OBB_System> Builder;
    Builder builder(systems,sel.max_nb_systems());
  
    T_Atom_classifier atom_classifier;
    std::vector<Vec3> points;

    Vec3 n_terminus_position;
    Vec3 c_terminus_position;

    std::vector<std::string> amino_acids{"ALA","ARG","ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS","ILE",
                                         "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};

    if (ESBTL::Line_reader<ESBTL::PDB::Line_format<>,Line_Selector,Builder>(sel,builder).template read_stream(std::cin,Accept_none_occupancy_policy())){
        if ( systems.empty() || systems[0].has_no_model() ){
            std::cerr << "No atoms found" << std::endl;
            return EXIT_FAILURE;
        }
        //Consider only the first model of the first system
        const OBB_System::Model& model=* systems[0].models_begin();

        for (OBB_System::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
            points.push_back(*it_atm);

            //store the N-terminus location
            if (it_atm->residue_sequence_number()==1 && it_atm->atom_name()=="CA"){

                double n_terminus_x = it_atm->x();
                double n_terminus_y = it_atm->y();
                double n_terminus_z = it_atm->z();

                n_terminus_position.Set(n_terminus_x, n_terminus_y, n_terminus_z);
            }

            //store the C -terminus location
            if (std::find(std::begin(amino_acids), std::end(amino_acids), it_atm->residue_name()) != std::end(amino_acids) && it_atm->atom_name()=="CA"){
                
                double c_terminus_x = it_atm->x();
                double c_terminus_y = it_atm->y();
                double c_terminus_z = it_atm->z();

                c_terminus_position.Set(c_terminus_x, c_terminus_y, c_terminus_z);
            }   
        }
    }
    else
        return EXIT_FAILURE;

    // build the minimum bounding box
    OBB obb_hull;
    obb_hull.build_from_convex_hull(points);

    // determine direction of rotation described in quaternion 
    //verifyRotation(obb_hull, points);

    // determine termini offsets from bounding box center
    Vec3 c_terminus_offset_rotated;
    c_terminus_offset_rotated = c_terminus_position - obb_hull.position();

    Vec3 n_terminus_offset_rotated;
    n_terminus_offset_rotated = n_terminus_position - obb_hull.position();

    // apply inverse rotation to the termini offsets to axis align
    Vec3 c_terminus_offset_aligned;
    c_terminus_offset_aligned = obb_hull.rotation().inverse()._transformVector(c_terminus_offset_rotated);

    Vec3 n_terminus_offset_aligned;
    n_terminus_offset_aligned = obb_hull.rotation().inverse()._transformVector(n_terminus_offset_rotated);

    // calculate the new positions
    Vec3 c_terminus_position_aligned;
    c_terminus_position_aligned = obb_hull.position() + c_terminus_offset_aligned;

    Vec3 n_terminus_position_aligned;
    n_terminus_position_aligned = obb_hull.position() + n_terminus_offset_aligned;

    //std::cout<<"Position: "<<obb_hull.position().transpose()<<std::endl;
    //std::cout<<"{\"extents\": ["<<obb_hull.extents()[0]<<", "<<obb_hull.extents()[1]<<", "<<obb_hull.extents()[2]<<"]}"<<std::endl;
    //std::cout<<"Rotation: "<<obb_hull.rotation().w()<<std::endl;

    // take the longest bounding box extent as the capsule's collision extent
    double capsule_collision_extent = std::max(std::max(obb_hull.extents()[0],obb_hull.extents()[1]),obb_hull.extents()[2]);

    // use the average of the remaining extents for the capsule radius
    double capsule_radius = (obb_hull.extents()[0] + obb_hull.extents()[1] + obb_hull.extents()[2] - capsule_collision_extent) * 0.5;

    std::cout<<capsule_radius<<std::endl;
    std::cout<<capsule_collision_extent<<std::endl;
    std::cout<<n_terminus_offset_aligned.transpose()<<std::endl;
    std::cout<<c_terminus_offset_aligned.transpose()<<std::endl;
    std::cout<<obb_hull.position().transpose()<<std::endl;


    return EXIT_SUCCESS;
}


void verifyRotation(OBB obb_hull, std::vector<Vec3> points){

    double center_x = obb_hull.position()[0];
    double center_y = obb_hull.position()[1];
    double center_z = obb_hull.position()[2];

    double extent_x = obb_hull.extents()[0];
    double extent_y = obb_hull.extents()[1];
    double extent_z = obb_hull.extents()[2];

    Vec3 points_forward;
    Vec3 points_inverse;

    //apply rotation and inverse rotation to all points
    for(int i = 0; i < points.size(); i++){

        points_forward = obb_hull.rotation()._transformVector(points[i]);

        // determine if it falls within the axis-aligned bounding box
        if((points_forward[0] <= center_x - extent_x*0.5 || points_forward[0] >= center_x + extent_x*0.5) && 
           (points_forward[1] <= center_y - extent_y*0.5 || points_forward[1] >= center_y + extent_y*0.5) &&
           (points_forward[2] <= center_z - extent_z*0.5 || points_forward[2] >= center_z + extent_z*0.5)) {

            std::cout<< "forward rotated point falls outside aligned box" << std::endl;
        }

        points_inverse = obb_hull.rotation().inverse()._transformVector(points[i]);

        if((points_inverse[0] <= center_x - extent_x*0.5 || points_inverse[0] >= center_x + extent_x*0.5) && 
           (points_inverse[1] <= center_y - extent_y*0.5 || points_inverse[1] >= center_y + extent_y*0.5) &&
           (points_inverse[2] <= center_z - extent_z*0.5 || points_inverse[2] >= center_z + extent_z*0.5)) {

            std::cout<< "inverse rotated point falls outside aligned box" << std::endl;
        }

    }
 
}
