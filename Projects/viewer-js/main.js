'use strict';

const electron=require('electron');
const app=electron.app;
const BrowserWindow=electron.BrowserWindow;
var ipc=electron.ipcMain;

var mainWindow=null;

app.on('window-all-closed', function(){
    if(process.platform!='darwin'){
        app.quit();
    }
});

// read scene file
var parsedData;
var inputChunks=[];
process.stdin.resume();
process.stdin.setEncoding('utf8');
process.stdin.on('data',function(chunk){
    inputChunks.push(chunk);
});
process.stdin.on('end',function(){
    var inputJSON=inputChunks.join();
    parsedData=JSON.parse(inputJSON);
});

ipc.on('synchronous-message',function(event,arg){
    event.returnValue=parsedData;
});


app.on('ready', function(){
    mainWindow=new BrowserWindow({width: 800, height: 600});
    mainWindow.loadURL('file://'+__dirname+'/index.html');

    mainWindow.on('closed',function(){
        mainWindow=null;
    });
});

