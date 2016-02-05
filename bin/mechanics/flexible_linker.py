from random_helpers import *
import math
import string

# handle flexible linkers
def flexible_linker(node,structure_nodes,force_nodes):
    links=node['links']
    subnodes=[]
    position=[x+y+z for x,y,z in zip(structure_nodes[node['structure1']]['position'],node['offset1'],[0,0,0])]
    constraints=[]
    distance=math.sqrt(sum([i**2 for i in node['link_offset']]))
    collision_extent=0
    subnode_offset=0
    if 'collision_extent' in node:
        collision_extent=node['collision_extent']
        subnode_offset=collision_extent+node['link_radius']

    for i in range(0,links):
        subnode={}
        position=[x+y+z for x,y,z in zip(position,node['link_offset'],[0,0,collision_extent])]
        subnode['type']='RIGID_STRUCTURE'
        subnode['radius']=node['link_radius']
        subnode['position']=position
        position=[x+y for x,y in zip(position,[0,0,collision_extent])]
        subnode['name']=random_name()
        subnode['collision_extent']=collision_extent
        subnodes.append(subnode)
        #structure_nodes.append(subnode)
        structure_nodes[subnode['name']]=subnode
        #expanded_data['root'].append(subnode)
    constraints.append({'structure1': node['structure1'],
                        'offset1': node['offset1'],
                        'structure2': subnodes[0]['name'],
                        'offset2': [0,0,-subnode_offset],
                        'distance': distance
                    })
    constraints.append({'structure1': subnodes[-1]['name'],
                        'offset1': [0,0,subnode_offset],
                        'structure2': node['structure2'],
                        'offset2': node['offset2'],
                        'distance': distance
                    })
    structure_nodes[node['structure2']]['position']=[x+y+z-t for x,y,z,t in zip(position,[0,0,collision_extent],node['link_offset'],node['offset2'])]
    for i in range(1,links):
        constraints.append({'structure1': subnodes[i-1]['name'],
                            'offset1': [0,0,subnode_offset],
                            'structure2': subnodes[i]['name'],
                            'offset2': [0,0,-subnode_offset],
                            'distance': distance})
    constraint={'type': 'RELATIVE_POSITION_CONSTRAINT','constraints': constraints}
    force_nodes.append(constraint)
