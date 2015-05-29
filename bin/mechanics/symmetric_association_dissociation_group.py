import math
import string
import json

def symmetric_association_dissociation_group(node,structure_nodes,force_nodes):
    interaction_types=node['interactions']
    interactions=[]
    for interaction_type in interaction_types:
        interaction={'bond_distance_threshold': interaction_type['bond_distance_threshold'],
                     'bond_orientation_threshold': interaction_type['bond_orientation_threshold'],
                     'base_dissociation_time': interaction_type['base_dissociation_time'],
                     'base_association_time': interaction_type['base_association_time'],
                     'binder_orientation': interaction_type['binder_orientation'],
                     'site_offset': interaction_type['site_offset']
                 }
        site_tag=interaction_type['site']['tag']
        site_offset=interaction_type['site']['site']
        interaction['sites']=[{'name': x['name'],'site': site_offset} for x in filter(lambda x: 'tag' in x and x['tag']==site_tag,structure_nodes.values())]
        interactions.append(interaction)
    constraint={'type': 'SYMMETRIC_ASSOCIATION_DISSOCIATION_CONSTRAINT', 'interactions': interactions}
    force_nodes.append(constraint)
