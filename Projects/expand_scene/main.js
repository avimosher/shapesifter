expand=require('./js/expand_scene');

var parsedData;
var inputChunks=[];
process.stdin.resume();
process.stdin.on('data',chunk => inputChunks.push(chunk));
process.stdin.on('end',function(){
    var inputJSON=JSON.parse(inputChunks.join());
    expand.expand(inputJSON,function(result){
        console.log(JSON.stringify(result,null,2));
    });
    //console.log(JSON.stringify(expand.expand(inputJSON),null,2));
});
