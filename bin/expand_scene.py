#!/usr/bin/python
import json
import sys
import random
import string
import math

json_data=open(sys.argv[1])
data=json.load(json_data)

expanded_data={}
expanded_data['root']=[]
expanded_data['output_directory']=data['output_directory']

def random_name():
    return ''.join(random.choice(string.ascii_lowercase+string.digits) for i in range(12))

for node in data['root']:
    if node['type']=='FLEXIBLE_LINKER':
        links=node['links']
        subnodes=[]
        position=[0,0,0]
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
            expanded_data['root'].append(subnode)
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
        for i in range(1,links):
            constraints.append({'structure1': subnodes[i-1]['name'],
                                'offset1': [0,0,0],
                                'structure2': subnodes[i]['name'],
                                'offset2': [0,0,0],
                                'distance': distance})
        constraint={'type': 'RELATIVE_POSITION_CONSTRAINT','constraints': constraints}
        expanded_data['root'].append(constraint)
    else:
        expanded_data['root'].append(node)

with open(sys.argv[2],'w') as outfile:
    json.dump(expanded_data,outfile,indent=4)

