if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['three'],function(THREE){
    var module={};

    module.create=function(node,structures,forces){
        var material = new THREE.LineBasicMaterial({color: 0xff0000, linewidth: 1});
        // loop over all constraints
        for (var j=0; j < node.constraints.length; j++){
            // retrieve structures and offsets
            s1 = node.constraints[j].structure1;
            s2 = node.constraints[j].structure2;
            n1 = new THREE.Vector3(node.constraints[j].offset1[0], node.constraints[j].offset1[1], node.constraints[j].offset1[2]);
            n2 = new THREE.Vector3(node.constraints[j].offset2[0], node.constraints[j].offset2[1], node.constraints[j].offset2[2]);

            first_attachment = new THREE.Vector3();
            first_attachment.addVectors(structures.getObjectByName(s1).position, n1);               

            second_attachment = new THREE.Vector3();
            second_attachment.addVectors(structures.getObjectByName(s2).position, n2);

            var lineGeometry = new THREE.Geometry();
            lineGeometry.vertices.push(first_attachment,second_attachment);

            var line = new THREE.Line(lineGeometry, material);
            line.name = [s1,s2,n1,n2];
            line.type="RELATIVE_POSITION_CONSTRAINT";
            forces.add(line);
        }
    };

    module.update=function(force,structures,forces){
        // retrieve relevant structures and offsets
        s1 = force.name[0];
        s2 = force.name[1];
        n1_copy = force.name[2].clone();
        n2_copy = force.name[3].clone();

        // rotate offset based on structure's current orientation   
        n1_copy.applyEuler(structures.getObjectByName(s1).rotation)
        n2_copy.applyEuler(structures.getObjectByName(s2).rotation)

        // add rotated offset to structure's center of mass
        force.geometry.vertices[0].addVectors(structures.getObjectByName(s1).position, n1_copy);
        force.geometry.vertices[1].addVectors(structures.getObjectByName(s2).position, n2_copy);
        force.geometry.verticesNeedUpdate = true;
    };
    return module;
});
