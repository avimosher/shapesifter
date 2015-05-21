import math
import string
import json

def association_dissociation_constraint(node,structure_nodes,force_nodes):
    interaction_types=node['interactions']
    interactions=[]
    for interaction_type in interaction_types:
        interaction={'bond_distance_threshold': interaction_type['bond_distance_threshold'],
                     'bond_orientation_threshold': interaction_type['bond_orientation_threshold'],
                     'base_dissociation_time': interaction_type['base_dissociation_time'],
                     'base_association_time': interaction_type['base_association_time'],
                     'relative_orientation': interaction_type['relative_orientation']}
        first_site_tag=interaction_type['first_site']['tag']
        first_site_offset=interaction_type['first_site']['site']
        second_site_tag=interaction_type['second_site']['tag']
        second_site_offset=interaction_type['second_site']['site']
        interaction['first_sites']=[{'name': x['name'],'site': first_site_offset} for x in filter(lambda x: 'tag' in x and x['tag']==first_site_tag,structure_nodes.values())]
        interaction['second_sites']=[{'name': x['name'],'site': second_site_offset} for x in filter(lambda x: 'tag' in x and x['tag']==second_site_tag,structure_nodes.values())]
        interactions.append(interaction)
    constraint={'type': 'ASSOCIATION_DISSOCIATION_CONSTRAINT', 'interactions': interactions}
    force_nodes.append(constraint)
