if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['three'],function(THREE){
    var module={};

    var orientation_quat = new THREE.Quaternion();
    var randomColor=require('randomcolor');

    function chooseColor(s){
        return randomColor(s);
    }

    module.create=function(node,data,structures,forces){
        // extract radius and capsule extent
        var radius = node.radius; 
        var capsule_extent = node.capsule_extent;

        var geometry;
        var material;

        // capsules
        if (capsule_extent > 0){
            // container object for capsule elements
            geometry = new THREE.Geometry();

            //var proteinMaterial = new THREE.MeshBasicMaterial( { map: texture, wireframe: false } );
            material = new THREE.MeshLambertMaterial( { wireframe: true, color: chooseColor(node.name) } );

            // define cylinder and two spheres to create a capsule
            var cylinder = new THREE.CylinderGeometry(radius, radius, capsule_extent*2, 17);
            var matrix = new THREE.Matrix4();
            matrix.makeRotationX(Math.PI/2);
            cylinder.applyMatrix(matrix);
            var top = new THREE.SphereGeometry(radius,17,17);
            var bottom = new THREE.SphereGeometry(radius,17,17);

            // set position of bottom and top spheres
            var m1 = new THREE.Matrix4();
            m1.makeTranslation(0,0,capsule_extent);
            top.applyMatrix(m1);

            var m2 = new THREE.Matrix4();
            m2.makeTranslation(0,0,-capsule_extent);
            bottom.applyMatrix(m2);

            // add spheres and cylinder to merged object
            geometry.merge(cylinder);
            geometry.merge(bottom);
            geometry.merge(top);
        }
        // spheres
        else{
            //var linkerMaterial = new THREE.MeshLambertMaterial( {color: 0xCCFF33, wireframe: false});
            material = new THREE.MeshLambertMaterial( {color: chooseColor(node.name), wireframe: true});
            geometry = new THREE.SphereGeometry(radius, 17, 17);
        }
        var structure = new THREE.Mesh(geometry, material);
        structure.name = node.name;
        structure.type = "RIGID_STRUCTURE";
        structures.add(structure);
    };

    module.update = function(frame_handle, structure){
        // retrieve position
        var updated_position = frame_handle.frame.position.data;

        // retrieve orientation
        var updated_orientation = frame_handle.frame.orientation.value0.data;

        // update position
        structure.position.x = updated_position[0];
        structure.position.y = updated_position[1];
        structure.position.z = updated_position[2];

        // set updated orientation
        orientation_quat.set(updated_orientation[0], updated_orientation[1], updated_orientation[2], updated_orientation[3]);

        // apply orientation to structure
        structure.rotation.setFromQuaternion(orientation_quat.normalize());
    };

    return module;
});
