if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['mathjs','three','./random_helpers'],function(math,THREE,random_helpers){
    return {handle: function(node,structureNodes,forceNodes){
        var links=node['links'];
        var color=node['color'];
        var subnodes=[];
        var structure1=structureNodes[node['structure1']];
        var structure2=structureNodes[node['structure2']];

        // correctly handle rotations with offsets
        var firstStructureOffset=node['offset1'].slice();
        if(orientation=structure1['orientation']){
            var axis3=new THREE.Vector3();axis3.fromArray(orientation.axis);
            var offset3=new THREE.Vector3();offset3.fromArray(firstStructureOffset);
            offset3.applyAxisAngle(axis3,orientation.angle*Math.PI/180);
            offset3.toArray(firstStructureOffset);
        }
        var secondStructureOffset=node['offset2'].slice();
        if(orientation=structure2['orientation']){
            var axis3=new THREE.Vector3();axis3.fromArray(orientation.axis);
            var offset3=new THREE.Vector3();offset3.fromArray(secondStructureOffset);
            offset3.applyAxisAngle(axis3,orientation.angle*Math.PI/180);
            offset3.toArray(secondStructureOffset);
        }
        var position=math.add(structure1['position'],firstStructureOffset);

        var constraints=[];
        var distance=math.sqrt(node['link_offset'].reduce(function(x,s){return x+s*s;},0));
        var collisionExtent=0;
        var subnodeOffset=0;
        if(node['collision_extent']){
            collisionExtent=node['collision_extent'];
            subnodeOffset=collisionExtent+node['link_radius'];
        }

        for(var i=0;i<links;i++){
            var subnode={};
            position=math.add(position,math.add(node['link_offset'],[0,0,collisionExtent]));
            //console.log(position);
            subnode['type']='RIGID_STRUCTURE';
            subnode['radius']=node['link_radius'];
            subnode['position']=position.slice();
            if(color){subnode['color']=color;}
            position=math.add(position,[0,0,collisionExtent]);
            subnode['name']=random_helpers.name();
            subnode['collision_extent']=collisionExtent;
            subnodes.push(subnode);
            structureNodes[subnode['name']]=subnode;
        }


        constraints.push({'structure1': node['structure1'],
                          'offset1': node['offset1'],
                          'structure2': subnodes[0]['name'],
                          'offset2': [0,0,-subnodeOffset],
                          'distance': distance});
        constraints.push({'structure1': subnodes[subnodes.length-1]['name'],
                          'offset1': [0,0,subnodeOffset],
                          'structure2': node['structure2'],
                          'offset2': node['offset2'],
                          'distance': distance});
        structureNodes[node['structure2']]['position']=math.add(position,math.add([0,0,collisionExtent],math.subtract(node['link_offset'],secondStructureOffset)));
        for(var i=1;i<links;i++){
            constraints.push({'structure1': subnodes[i-1]['name'],
                              'offset1': [0,0,subnodeOffset],
                              'structure2': subnodes[i]['name'],
                              'offset2': [0,0,-subnodeOffset],
                              'distance': distance});}
        forceNodes.push({'type': 'RELATIVE_POSITION_CONSTRAINT','constraints': constraints});
    }}});
