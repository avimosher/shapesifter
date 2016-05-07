if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(function(){

    var module={};

    module.expand=function(data){
        var expandedData={};
        expandedData['root']=[];
        for(key in data){
            if(key!='root'){
                expandedData[key]=data[key];
            }
        }
        var structureNodes={};
        var forceNodes=[];

        function rigidStructure(node,structureNodes,forceNodes){
            structureNodes[node['name']]=node;}
        function force(node,structureNodes,forceNodes){
            forceNodes.push(node);}
        
        function handle(node,structureNodes,forceNodes){
            try{
                require('./'+node['type']).handle(node,structureNodes,forceNodes);
            }catch(err){
                //console.log(err);
                //console.log('could not find '+node['type']);
                expandedData['root'].push(node);}}

        data['root'].forEach(function(node){
            handle(node,structureNodes,forceNodes);
        });

        for(key in structureNodes){
            expandedData['root'].push(structureNodes[key]);}
        Array.prototype.push.apply(expandedData['root'],forceNodes);
        return expandedData;
    };
    return module;
});