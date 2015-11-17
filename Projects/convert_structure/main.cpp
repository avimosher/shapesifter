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

char* Get_Command_Option(char** begin,char** end,const std::string& option)
{
    char **itr=std::find(begin,end,option);
    if(itr!=end){
        if(++itr!=end){return *itr;}
        else{return *(--itr);}}
    return 0;
}

int main(int argc,char **argv)
{
    typedef double T;
    typedef Matrix<T,3,1> TV;
    typedef ESBTL::Molecular_system<ESBTL::Default_system_items,Vec3> OBB_System;
    auto simulation=std::make_shared<SIMULATION<TV>>();

    char *restart=Get_Command_Option(argv,argv+argc,"-restart");
    if(restart){simulation->Set_Restart(std::stoi(restart));}
    if(Get_Command_Option(argv,argv+argc,"-substeps")){
        simulation->substeps=true;
    }
    if(Get_Command_Option(argv,argv+argc,"-log")){
        LOG::Output(true);
    }
    if(Get_Command_Option(argv,argv+argc,"-nowrite")){
        simulation->write=false;
    }



    if (argc != 2 ){
        std::cerr << "Please provide a filename"  << std::endl;  
        return EXIT_FAILURE;
    }
  
    ESBTL::PDB_line_selector_two_systems sel;
  
    std::vector<OBB_System> systems;
    ESBTL::All_atom_system_builder<OBB_System> builder(systems,sel.max_nb_systems());
  
    T_Atom_classifier atom_classifier;
    std::vector<Vec3> points;

    if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
        if ( systems.empty() || systems[0].has_no_model() ){
            std::cerr << "No atoms found" << std::endl;
            return EXIT_FAILURE;
        }
        //Consider only the first model of the first system
        const OBB_System::Model& model=* systems[0].models_begin();
        for (OBB_System::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
            points.push_back(*it_atm);
            //std::cout << it_atm->x() << " " << it_atm->y() << " " << it_atm->z() << " " ;
            //std::cout << atom_classifier.get_properties(*it_atm).value() << "\n";
        }
    }
    else
        return EXIT_FAILURE;

    OBB obb_hull;
    obb_hull.build_from_convex_hull(points);
    std::cout<<"OBB convex hull volume: "<<obb_hull.volume()<<std::endl;
    std::cout<<"Position: "<<obb_hull.position().transpose()<<std::endl;
    std::cout<<"Extents: "<<obb_hull.extents().transpose()<<std::endl;
    std::cout<<"Rotation: "<<obb_hull.rotation()<<std::endl;
    return EXIT_SUCCESS;
}
