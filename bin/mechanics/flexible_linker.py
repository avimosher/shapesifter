import random
import math
import string

def random_name():
    return ''.join(random.choice(string.ascii_lowercase+string.digits) for i in range(12))
def random_vector(low,high):
    return [random.uniform(l,h) for l,h in zip(low,high)]

# handle flexible linkers
def flexible_linker(node,structure_nodes,force_nodes):
    links=node['links']
    subnodes=[]
    position=[x+y+z for x,y,z in zip(structure_nodes[node['structure1']]['position'],node['offset1'],[0,0,0])]
    constraints=[]
    distance=math.sqrt(sum([i**2 for i in node['link_offset']]))
    for i in range(0,links):
        subnode={}
        position=[x+y for x,y in zip(position,node['link_offset'])]
        subnode['type']='RIGID_STRUCTURE'
        subnode['radius']=node['link_radius']
        subnode['position']=position
        subnode['name']=random_name()
        subnodes.append(subnode)
        #structure_nodes.append(subnode)
        structure_nodes[subnode['name']]=subnode
        #expanded_data['root'].append(subnode)
    constraints.append({'structure1': node['structure1'],
                        'offset1': node['offset1'],
                        'structure2': subnodes[0]['name'],
                        'offset2': [0,0,0],
                        'distance': distance
                    })
    constraints.append({'structure1': subnodes[-1]['name'],
                        'offset1': [0,0,0],
                        'structure2': node['structure2'],
                        'offset2': node['offset2'],
                        'distance': distance
                    })
    structure_nodes[node['structure2']]['position']=[x+y-z for x,y,z in zip(subnodes[-1]['position'],node['link_offset'],node['offset2'])]
    for i in range(1,links):
        constraints.append({'structure1': subnodes[i-1]['name'],
                            'offset1': [0,0,0],
                            'structure2': subnodes[i]['name'],
                            'offset2': [0,0,0],
                            'distance': distance})
    constraint={'type': 'RELATIVE_POSITION_CONSTRAINT','constraints': constraints}
    force_nodes.append(constraint)
