if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['three'],function(THREE){

    var initScene;
    var frame_directory='';
    var scene;
    var camera;
    var controls;
    var sceneInitialized=false;
    var structures=new THREE.Object3D();
    var forces=new THREE.Object3D();

    // prep animation
    var frame_index = 0;
    var running;
    var writing;
    var titles = true;
    var frame_total = 0;


    function readSceneData(callback, scene){
        const fs=require('fs');
        try{
            fs.openSync(scene,'r+');
            callback(JSON.parse(fs.readFileSync(scene)));}
        catch(err){
            console.log('Failed to read' + scene);}
    }

    function createStructuresForces(scene_data, structures, forces) {
        for (var i = 0; i < scene_data.root.length; i++){
            try{require('./'+scene_data.root[i].type).create(scene_data.root[i],structures,forces);}
            catch(e){
                //console.log(scene_data.root[i].type);
                //console.log(e);
            }}
        /*var loader=new THREE.OBJLoader();
          loader.load('untitled.obj',function(object){
          structures.add(object);
          });*/
    }

    function updateStructures(index, structures){
        const fs=require('fs');
        try{
            var frame_data=JSON.parse(fs.readFileSync(frame_directory+'/frame.'+index));
            for (var i = 0; i < structures.children.length; i++){
                frame_handle = frame_data.value2[0].ptr_wrapper.data.value0[i].ptr_wrapper.data;
                var structure = structures.getObjectByName(frame_handle.name);
                try{
                    require('./'+structure.type).update(frame_handle, structure);}
                catch(e){console.log(e);}}}
        catch(err){console.log(err);}
    }

    function updateForces(forces, structures){
        for (var i = 0; i < forces.children.length; i++){
            try{ require('./'+forces.children[i].type).update(forces.children[i],structures,forces);}
            catch(e){console.log(e);}}}

    function countFrames(callback){
        // find number of frames in frame folder
        var fs=require('fs');
        callback(fs.readdirSync(frame_directory).length);
    }

    var initializeViewer = function() {
        // set the scene size
        var WIDTH = window.innerWidth,HEIGHT = window.innerHeight;

        // set some camera attributes
        var VIEW_ANGLE = 45,
                ASPECT = WIDTH / HEIGHT,
                NEAR = 0.1,
                FAR = 10000;

        // create a WebGL renderer, camera, and a scene
        var canvas = document.getElementById("protein_viewer");
        renderer = new THREE.WebGLRenderer({ canvas: canvas });

        camera = new THREE.PerspectiveCamera( VIEW_ANGLE, WIDTH / HEIGHT, NEAR, FAR );
        scene = new THREE.Scene();

        var bounds = new THREE.Box3();
        bounds.setFromObject(structures);
        camera.position.addVectors(bounds.center(),new THREE.Vector3(0,0,1.5*bounds.size().length()));

        var OrbitControls=require('three-orbit-controls')(THREE);
        controls = new OrbitControls(camera, renderer.domElement);
        controls.target.copy(bounds.center());
        controls.update();

        // add the camera to the scene
        scene.add(camera);

        // start the renderer
        renderer.setSize(WIDTH, HEIGHT);
        // document.body.append( renderer.domElement );

        // lights
        var forwardDirectionalLight = new THREE.DirectionalLight(0xffffff);
        forwardDirectionalLight.position.set(1,1,1);
        scene.add(forwardDirectionalLight);

        var reverseDirectionalLight = new THREE.DirectionalLight(0xffffff);
        reverseDirectionalLight.position.set(-1,-1,-1);
        scene.add(reverseDirectionalLight);

        var ambientLight = new THREE.AmbientLight(0x222222);
        scene.add(ambientLight);

        // read scene data
        // var scene_data;
        /*readSceneData(function (parsed_scene){
          scene_data = parsed_scene;
          }, scene_file);*/
        window.onkeydown=function(e){
            var continueProcessing=false;
            switch(e.keyCode){
                case 189: // - - back
                incrementFrame(-1);
                break;

                case 187: // + - forward
                incrementFrame(1);
                break;

                case 27: // esc - quit
                var ipc=require('electron').ipcRenderer;
                ipc.sendSync('synchronous-message','quit');
                break;

                case 71: // g - go to frame
                var dialogPolyfill=require('dialog-polyfill');
                var dialog=$('<dialog>');
                dialog.append($('<p>',{'text': 'Go to frame'}));
                var text=$('<input>',{'type': 'text'});
                var goToFrame=function(){
                    dialog[0].close();
                    console.log(text.val());
                    frame_index=parseInt(text.val());
                    incrementFrame(0);
                    dialog[0].remove();
                };
                text.keypress(function(e){
                        if(e.which==13){
                            goToFrame();
                            return false;}});
                dialog.append(text);
                var button=$('<button>',{'text': 'Go'});
                button.on('click',goToFrame);
                dialog.append($('<p>').append(button));
                $('body').append(dialog);
                dialog[0].showModal();
                //$("#dialog")[0].showModal();
                
                /*dialog.open({message: 'Go to frame',
                             input: "<input name=\"frame\" type=\"text\" placeholder=\"\" required />\n",
                             callback: function(frame){
                                 alert(frame);
                                 frame_index=frame;
                                 incrementFrame(0);
                                 }});*/
                break;

                case 80: // p - toggle play
                running=!running;
                break;

                case 82: // r - reset
                running=false;
                frame_index=0;
                incrementFrame(0);
                break;

                case 84: // t - titles
                titles=!titles;
                if(titles){
                    window.$("#label").show();}
                else{
                    window.$("#label").hide();}
                break;

                case 87: // w - toggle writing
                writing=!writing;
                break;
                
                default:
                    continueProcessing=true;
                break;
            };
            return continueProcessing;
        };
    }

    function update() {
        // update structures
        updateStructures(frame_index, structures);

        // update forces
        updateForces(forces, structures);
    }

    function incrementFrame(increment){
        frame_index += increment;
        if (frame_index >= frame_total){
            frame_index = frame_total - 1;}
        if (frame_index < 0){frame_index = 0;}
        window.$("#label").text("Frame "+frame_index);
        update();
    }

    function render() {
        if (running) {incrementFrame(1);}

        renderer.render(scene, camera);
        controls.update();
        if(!writing){requestAnimationFrame(render);}
        else{
            try{
                var remote=require('remote');
                remote.getCurrentWindow().capturePage(function(buf){
                    var S=require('string');
                    remote.require('fs').writeFile(frame_directory+'/screenshot.'+S(frame_index).padLeft(5,'0')+'.png',buf.toPng(),function(){
                        requestAnimationFrame(render);
                    })});}
            catch(err){}}
    };

    var module={};
    module.main = function() {
        if(sceneInitialized){
            scene.remove(structures);
            scene.remove(forces);}
        else{
            var ipc=require('electron').ipcRenderer;
            scene_data=ipc.sendSync('synchronous-message','');}

        frame_directory=scene_data['output_directory'];

        // create structures and forces
        createStructuresForces(scene_data, structures, forces);

        // count total frames
        frame_total = 0;
        countFrames(function (frames){
            frame_total += frames;
        });

        // prep animation
        frame_index = 0;
        running = false;
        writing = false;

        incrementFrame(0);

        if(!sceneInitialized){
            initializeViewer();
            sceneInitialized=true;
        }
        // add the structures and forces
        scene.add(structures);
        scene.add(forces);

        render();
    };
    return module;
});

//window.onload = initScene;
//window.addEventListener( 'load', main, true );
