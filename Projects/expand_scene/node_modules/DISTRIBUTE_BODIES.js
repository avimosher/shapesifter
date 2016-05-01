if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['random_helpers'],function(random_helpers){
    return {handle: function(node,structureNodes,forceNodes){
        bodies=node['bodies'];
        for(var i=0;i<bodies;i++){
            subnode={};
            subnode['type']='RIGID_STRUCTURE';
            subnode['radius']=node['radius'];
            subnode['position']=random_helpers.vector(node['min'],node['max']);
            subnode['name']=random_helpers.name();
            if(node['tag']){
                subnode['tag']=node['tag'];
            }
            structureNodes[subnode['name']]=subnode;
        }
    }
           }});
