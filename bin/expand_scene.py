#!/usr/bin/python
import json
import sys
import string
from mechanics.flexible_linker import *
from mechanics.association_dissociation_group import *
from mechanics.symmetric_association_dissociation_group import *
from mechanics.distribute_bodies import *

data=json.load(sys.stdin)

expanded_data={}
expanded_data['root']=[]
for key in data:
    if key!='root':
        expanded_data[key]=data[key]


structure_nodes={}
force_nodes=[]

def rigid_structure(node,structure_nodes,force_nodes):
    structure_nodes[node['name']]=node
def force(node,structure_nodes,force_nodes):
    force_nodes.append(node)


handle_node={'FLEXIBLE_LINKER': flexible_linker,
             'DISTRIBUTE_BODIES': distribute_bodies,
             'ASSOCIATION_DISSOCIATION_CONSTRAINT': force,
             'ASSOCIATION_DISSOCIATION_GROUP': association_dissociation_group,
             'SYMMETRIC_ASSOCIATION_DISSOCIATION_GROUP': symmetric_association_dissociation_group,
             'RIGID_STRUCTURE': rigid_structure,
             'RELATIVE_POSITION_CONSTRAINT': force,
             'ABSOLUTE_POSITION_CONSTRAINT': force}

def handle_default(node,structure_nodes,force_nodes):
    expanded_data['root'].append(node)


for node in data['root']:
    handle_node.get(node['type'],handle_default)(node,structure_nodes,force_nodes)

def strip_tag(node):
    node.pop('tag',None)
    return node

expanded_data['root'].extend([strip_tag(obj) for obj in structure_nodes.values()])
expanded_data['root'].extend(force_nodes)

json.dump(expanded_data,sys.stdout,indent=4)

