if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(function(){
    return {handle: function(node,structureNodes,forceNodes){
        structureNodes[node['name']]=node;}}});
