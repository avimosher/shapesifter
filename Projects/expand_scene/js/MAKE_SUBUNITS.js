if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['js/FLEXIBLE_LINKER'],function(flexible_linker){
    return {handle: function(node,structureNodes,forceNodes){
        pore_radius = 2;
	capsule_radius = 0.8;
	capsule_extent = 3;
	number = 7;
        axis_x=5;;
        axis_y=10;
        axis_z=3;
	
	vector_magnitude = Math.sqrt((axis_x*axis_x) + (axis_y*axis_y) + (axis_z*axis_z));
        axis_x = axis_x/vector_magnitude;
        axis_y = axis_y/vector_magnitude;
        axis_z = axis_z/vector_magnitude;
	theta = (360/(number)) * (Math.PI/180);

	for(var i=0;i<number;i++){
            ang = i*theta;
            x =(pore_radius)*Math.cos(ang);
            y =(pore_radius)*Math.sin(ang);
            structureNodes['capsule '+i] = {'type': 'RIGID_STRUCTURE', 'kinematic': 1, 'position': [x,y,0], 'orientation': {'axis': [axis_x,axis_y,axis_z], 'angle': 0}, 'radius': capsule_radius, 'capsule_extent': capsule_extent, 'name': 'capsule '+i};
        }

        switcher={'low': 1./2,'med': 1,'high': 2};
        multiplier1=switcher[node['linker_length'][0]];
        multiplier2=switcher[node['linker_length'][1]];
        structureNodes['SC'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(7.9*multiplier1)], 'orientation': {'axis': [1, 0, 0], 'angle': 0}, 'radius': 1.95, 'capsule_extent': 0, 'name': 'SC'};
        structureNodes['phi29'] = {'type': 'RIGID_STRUCTURE', 'position': [0, 0,(14.6*multiplier1)], 'radius': 2.25, 'capsule_extent': 3, 'name': 'phi29'};
        fnode1 = {"type": "FLEXIBLE_LINKER", "structure1": "capsule 0", "offset1": [-1,0,2.05], "structure2": "SC", "offset2": [0,0,-1.95], "links": 3, "link_radius": (0.125*multiplier2), "link_offset": [0,0,(0.4*multiplier2)]};
        fnode2 = {"type": "FLEXIBLE_LINKER", "structure1": "SC", "offset1": [0,0,1.95], "structure2": "phi29", "offset2": [0,2.25,0], "links": 4, "link_radius": (0.125*multiplier2), "link_offset": [0,0,(0.4*multiplier2)]};
        flexible_linker.handle(fnode1,structureNodes,forceNodes);
        flexible_linker.handle(fnode2,structureNodes,forceNodes);
    }};
});


