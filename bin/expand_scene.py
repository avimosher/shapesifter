#!/usr/bin/python
import json
import sys
from pprint import pprint

json_data=open(sys.argv[1])
data=json.load(json_data)

expanded_data={}
expanded_data['root']=[]
expanded_data['output_directory']=data['output_directory']

for node in data['root']:
    if node['type']=='FLEXIBLE_LINKER':
        
    else:
        expanded_data['root'].append(node)

with open(sys.argv[2],'w') as outfile:
    json.dump(expanded_data,outfile,indent=4)

