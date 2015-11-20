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
#include "Vec3.h"
#include "OBB.h"
using namespace Mechanics;

typedef ESBTL::Generic_classifier<ESBTL::Radius_of_atom<double,ESBTL::Default_system::Atom> >                                  T_Atom_classifier;
typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

int main(int argc,char **argv)
{
    typedef double T;
    typedef ESBTL::Molecular_system<ESBTL::Default_system_items,Vec3> OBB_System;
    
    typedef ESBTL::PDB_line_selector_two_systems Line_Selector;
    Line_Selector sel;
  
    std::vector<OBB_System> systems;
    typedef ESBTL::All_atom_system_builder<OBB_System> Builder;
    Builder builder(systems,sel.max_nb_systems());
  
    T_Atom_classifier atom_classifier;
    std::vector<Vec3> points;

    if (ESBTL::Line_reader<ESBTL::PDB::Line_format<>,Line_Selector,Builder>(sel,builder).template read_stream(std::cin,Accept_none_occupancy_policy())){
        if ( systems.empty() || systems[0].has_no_model() ){
            std::cerr << "No atoms found" << std::endl;
            return EXIT_FAILURE;
        }
        //Consider only the first model of the first system
        const OBB_System::Model& model=* systems[0].models_begin();
        for (OBB_System::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
            points.push_back(*it_atm);
        }
    }
    else
        return EXIT_FAILURE;

    OBB obb_hull;
    obb_hull.build_from_convex_hull(points);
    //std::cout<<"Position: "<<obb_hull.position().transpose()<<std::endl;
    std::cout<<"{\"extents\": ["<<obb_hull.extents()[0]<<", "<<obb_hull.extents()[1]<<", "<<obb_hull.extents()[2]<<"]}"<<std::endl;
    //std::cout<<"Rotation: "<<obb_hull.rotation()<<std::endl;
    return EXIT_SUCCESS;
}
