from random_helpers import *

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

