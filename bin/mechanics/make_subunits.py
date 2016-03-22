import re
import math
import json
from flexible_linker import *


def make_subunits(node,structure_nodes,force_nodes):
	pore_radius = 2
	capsule_radius = 0.8
	collision_extent = 3
	number = 7
	center = "[5,10,3]"
	axis = "[5,10,3]"

	center_x, center_y, center_z = float(re.split("\s+|,|\[|]", center)[1]), float(re.split("\s+|,|\[|]", center)[2]), float(re.split("\s+|,|\[|]", center)[3])
	axis_x, axis_y, axis_z = float(re.split("\s+|,|\[|]", axis)[1]), float(re.split("\s+|,|\[|]", axis)[2]), float(re.split("\s+|,|\[|]", axis)[3])
	
	vector_magnitude = math.sqrt((axis_x**2) + (axis_y**2) + (axis_z**2))
        axis_x = axis_x/vector_magnitude
        axis_y = axis_y/vector_magnitude
        axis_z = axis_z/vector_magnitude
	theta = (360/(number)) * (math.pi/180)

	for i in range(number):
                ang = i*theta
                x =(pore_radius)*math.cos(ang)
                y =(pore_radius)*math.sin(ang)
                structure_nodes['capsule '+str(i)] = {'type': 'RIGID_STRUCTURE', 'kinematic': 1, 'position': [x,y,0], 'orientation': {'axis': [axis_x,axis_y,axis_z], 'angle': 0}, 'radius': capsule_radius, 'collision_extent': collision_extent, 'name': 'capsule '+str(i)}


        switcher={'low': 1./2,'med': 1,'high': 2}
        multiplier1=switcher[node['linker_length'][0]]
        multiplier2=switcher[node['linker_length'][1]]
        structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(7.9*multiplier1)], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
        structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6*multiplier1)], 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
        #structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6*multiplier1)], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
        fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [-1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": (0.125*multiplier2), "link_offset": [0,0,(0.4*multiplier2)]}
        fnode2 = {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": (0.125*multiplier2), "link_offset": [0,0,(0.4*multiplier2)]}
        flexible_linker(fnode1,structure_nodes,force_nodes)
        flexible_linker(fnode2,structure_nodes,force_nodes)
