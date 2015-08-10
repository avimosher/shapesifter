import re
import math
import json
from mechanics.flexible_linker import *


def make_subunits(node,structure_nodes,force_nodes):
	pore_radius = 2
	capsule_radius = 0.95
	collision_extent = 3
	number = 7
	center = "[5,10,3]"
	angle = 0
	axis = "[5,10,3]"

	center_x, center_y, center_z = float(re.split("\s+|,|\[|]", center)[1]), float(re.split("\s+|,|\[|]", center)[2]), float(re.split("\s+|,|\[|]", center)[3])

	axis_x, axis_y, axis_z = float(re.split("\s+|,|\[|]", axis)[1]), float(re.split("\s+|,|\[|]", axis)[2]), float(re.split("\s+|,|\[|]", axis)[3])
	
	vector_magnitude = math.sqrt((axis_x**2) + (axis_y**2) + (axis_z**2))
	if math.sqrt((axis_x**2) + (axis_y**2) + (axis_z**2)) != 1:
		axis_x = axis_x/vector_magnitude
		axis_y = axis_y/vector_magnitude
		axis_z = axis_z/vector_magnitude

	theta = (360/(number)) * (math.pi/180)
	rho = math.sqrt((center_x**2) + (center_y**2) + (center_z**2))
	phi = math.acos((axis_z**2)/(math.sqrt((axis_x**2) + (axis_y**2)))) * (180/math.pi)

	#print unit_x, unit_y, unit_z
	capsule_constraints =[]
	linker_constraints = []
	x_array = []
	y_array = []
	z_array = []
	for i in range(number):
		theta_multiplier = 1
		ang = 0
		while ang < 360:
			
			#print x_array, y_array, z_array

			ang = theta_multiplier*theta

			x =(pore_radius)*math.cos(ang)
			y =(pore_radius)*math.sin(ang)

			x_array.append(x)
			y_array.append(y)
			z_array.append(0)
			theta_multiplier = theta_multiplier + 1
			if len(x_array) == number:
				for j in range(number):


					structure_nodes['capsule '+str(j)] = {'type': 'RIGID_STRUCTURE', 'position': [x_array[j],y_array[j],z_array[j]], 'orientation': {'axis': [axis_x,axis_y,axis_z], 'angle': 0}, 'radius': capsule_radius, 'collision_extent': collision_extent, 'name': 'capsule '+str(j)}

					capsule_constraints.append({'type': 'linear', 'structure': 'capsule '+str(j), 'direction': [1,0,0], 'magnitude': x_array[j]})
					capsule_constraints.append({'type': 'linear', 'structure': 'capsule '+str(j), 'direction': [0,1,0], 'magnitude': y_array[j]})
					capsule_constraints.append({'type': 'linear', 'structure': 'capsule '+str(j), 'direction': [0,0,1], 'magnitude': z_array[j]})
					capsule_constraints.append({'structure': 'capsule '+str(j), 'type': 'angular', 'orientation': {'angle': 0, 'axis': [1,0,0]}})

				capsule_constraint ={"type": "ABSOLUTE_POSITION_CONSTRAINT", "constraints": capsule_constraints}
				force_nodes.append(capsule_constraint)


	if node['linker_length'][0] == 'low' and node['linker_length'][1] == 'low':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(7.9/2)], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6/2)], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": (0.125/2), "link_offset": [0,0,(0.4/2)]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": (0.125/2), "link_offset": [0,0,(0.4/2)]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)


	if node['linker_length'][0] == 'low' and node['linker_length'][1] == 'med':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(7.9/2)], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6/2)], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": 0.125, "link_offset": [0,0,0.4]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": 0.125, "link_offset": [0,0,0.4]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)


	if node['linker_length'][0] == 'low' and node['linker_length'][1] == 'high':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(7.9/2)], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6/2)], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": (0.125*2), "link_offset": [0,0,(0.4*2)]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": (0.125*2), "link_offset": [0,0,(0.4*2)]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)


	if node['linker_length'][0] == 'med' and node['linker_length'][1] == 'low':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,7.9], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,14.6], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": (0.125/2), "link_offset": [0,0,(0.4/2)]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": (0.125/2), "link_offset": [0,0,(0.4/2)]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)


	if node['linker_length'][0] == 'med' and node['linker_length'][1] == 'med':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,7.9], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,14.6], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": 0.125, "link_offset": [0,0,0.4]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": 0.125, "link_offset": [0,0,0.4]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)


	if node['linker_length'][0] == 'med' and node['linker_length'][1] == 'high':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,7.9], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,14.6], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": (0.125*2), "link_offset": [0,0,(0.4*2)]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": (0.125*2), "link_offset": [0,0,(0.4*2)]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)


	if node['linker_length'][0] == 'high' and node['linker_length'][1] == 'low':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(7.9*2)], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6*2)], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": (0.125/2), "link_offset": [0,0,(0.4/2)]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": (0.125/2), "link_offset": [0,0,(0.4/2)]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)



	if node['linker_length'][0] == 'high' and node['linker_length'][1] == 'med':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(7.9*2)], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6*2)], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": 0.125, "link_offset": [0,0,0.4]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": 0.125, "link_offset": [0,0,0.4]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)


	if node['linker_length'][0] == 'high' and node['linker_length'][1] == 'high':
		structure_nodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(7.9*2)], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'collision_extent': 0, 'name': 'SC'}
		structure_nodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6*2)], 'orientation': {'axis': [1, 0, 0], 'angle': -90}, 'radius': 2.25, 'collision_extent': 3, 'name': 'phi29'}
		fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [1.3,0,1.95], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": (0.125*2), "link_offset": [0,0,(0.4*2)]}

		fnode2 =         {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": (0.125*2), "link_offset": [0,0,(0.4*2)]}
		flexible_linker(fnode1,structure_nodes,force_nodes)
		flexible_linker(fnode2,structure_nodes,force_nodes)
