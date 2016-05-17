if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['requirejs'],function(requirejs){

    var module={};

    module.expand=function(data, callback){
        var expandedData={};
        expandedData['root']=[];
        for(key in data){
            if(key!='root'){
                expandedData[key]=data[key];
            }
        }
        var structureNodes={};
        var forceNodes=[];

        function handle(node,structureNodes,forceNodes){
            try{
                requirejs('./js/'+node['type']).handle(node,structureNodes,forceNodes,expandedData['root']);
            }catch(err){
                //console.log(err);
                //console.log('could not find '+node['type']);
                expandedData['root'].push(node);}}

        var types={};
        data['root'].forEach(function(node){
            types['./js/'+node['type']]=1;
        });
        requirejs(Object.keys(types),function(){
            data['root'].forEach(function(node){
                handle(node,structureNodes,forceNodes);
            });
            for(key in structureNodes){
                expandedData['root'].push(structureNodes[key]);}
            Array.prototype.push.apply(expandedData['root'],forceNodes);
            callback(expandedData);
        });
    };
    return module;
});
