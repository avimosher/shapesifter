expand=require('expand_scene');

var parsedData;
var inputChunks=[];
process.stdin.resume();
process.stdin.on('data',chunk => inputChunks.push(chunk));
process.stdin.on('end',function(){
    var inputJSON=JSON.parse(inputChunks.join());
    console.log(JSON.stringify(expand.expand(inputJSON),null,2));
});
