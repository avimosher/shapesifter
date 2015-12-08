#!/usr/bin/python

import json
import sys

data={}
data['root']=[]

structure_nodes={}
force_nodes=[]

def rigid_structure_node(structure_nodes,force_nodes,name,position,radius,extent):
    structure_nodes[name]={'type': 'RIGID_STRUCTURE',
                           'radius': radius,
                           'collision_extent': extent,
                           'name': name,
                           'position': position}
    

rigid_structure_node(structure_nodes,force_nodes,'test',[1,1,1],1,0)

data['root'].extend(structure_nodes.values())
data['root'].extend(force_nodes)

json.dump(data,sys.stdout,indent=4)
