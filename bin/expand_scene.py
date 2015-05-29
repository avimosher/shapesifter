#!/usr/bin/python
import json
import sys
import random
import string
from mechanics.flexible_linker import *
from mechanics.association_dissociation_group import *
from mechanics.symmetric_association_dissociation_group import *

data=json.load(sys.stdin)

expanded_data={}
expanded_data['root']=[]
for key in data:
    if key!='root':
        expanded_data[key]=data[key]


def random_name():
    return ''.join(random.choice(string.ascii_lowercase+string.digits) for i in range(12))
def random_vector(low,high):
    return [random.uniform(l,h) for l,h in zip(low,high)]

structure_nodes={}
force_nodes=[]

def rigid_structure(node,structure_nodes,force_nodes):
    structure_nodes[node['name']]=node


def distribute_bodies(node,structure_nodes,force_nodes):
    bodies=node['bodies']
    for i in range(0,bodies):
        subnode={}
        subnode['type']='RIGID_STRUCTURE'
        subnode['radius']=node['radius']
        subnode['position']=random_vector(node['min'],node['max'])
        subnode['name']=random_name()
        if 'tag' in node:
            subnode['tag']=node['tag']
        structure_nodes[subnode['name']]=subnode

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

#with open(sys.argv[2],'w') as outfile:
#    json.dump(expanded_data,outfile,indent=4)
json.dump(expanded_data,sys.stdout,indent=4)

