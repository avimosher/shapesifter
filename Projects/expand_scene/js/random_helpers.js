if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['chance'],function(Chance){
    var module={};
    module.chance=new Chance(812321);

    module.name=function(){return module.chance.string();};
    module.vector=function(low,high){
        var v=[];
        for(var i=0;i<3;i++){
            v[i]=module.chance.floating({min: low[i],max: high[i]});}
        return v;}
    return module;
});
