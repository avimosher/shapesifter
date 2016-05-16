if (typeof define !== 'function'){
    var define=require('amdefine')(module);
}

define(['three'],function(THREE){
    var module={};

    function createLine(material,points){
        //var lineGeometry=new THREE.Geometry();
        //var beginPoint=new THREE.Vector3();
        //var endPoint=new THREE.Vector3();
        var extent=[0,0,0];
        var epsilon=.01;
        for(var i=0;i<3;i++){
            extent[i]=points[1][i]-points[0][i]+2*epsilon;
        }
        console.log(extent);
        var boxGeometry=new THREE.BoxGeometry(extent[0],extent[1],extent[2]);
        boxGeometry.translate(points[0][0]+extent[0]/2+epsilon,points[0][1]+extent[1]/2+epsilon,points[0][2]+extent[2]/2+epsilon);
        return new THREE.Mesh(boxGeometry,material);
        //lineGeometry.vertices.push(beginPoint,endPoint);
        //return new THREE.Line(lineGeometry,material);
    }

    // TODO: change pass-through to take data
    module.create=function(node,data,structures,forces){
        //var material=new THREE.LineBasicMaterial({color: 0xffffff, linewidth: 1});
        var material=new THREE.MeshPhongMaterial({color: 0xffffff, wireframe: false});
        console.log('sup');
        var group=new THREE.Object3D();
        var walls=[];
        walls[0]=node.lower;
        walls[1]=node.upper;
        group.type="WALL_CONSTRAINT";
        var horizontalPoints=[[0,0,0],[0,0,0]];
        var verticalPoints=[[0,0,0],[0,0,0]];
        var edgeLengths=[];
        var domain=[[0,0,0],[0,0,0]];
        var d=3;
        for (var axis=0;axis<d;axis++){
            domain[0][axis]=data.domain.minimum_corner[axis];
            domain[1][axis]=data.domain.maximum_corner[axis];
            edgeLengths[axis]=domain[1][axis]-domain[0][axis];}
        for (var axis=0;axis<d;axis++){
            for (var w=0,sgn=1;w<2;w++,sgn-=2){
                if(walls[w][axis]){
                    var segments=6;
                    for (var i=0;i<=segments;i++){
                        // walk across each dimension in fractions
                        for(var j=0;j<2;j++){
                            horizontalPoints[j][axis]=domain[w][axis];
                            horizontalPoints[j][(axis+1)%d]=domain[0][(axis+1)%d]+i*1.0/segments*edgeLengths[(axis+1)%d];
                            horizontalPoints[j][(axis+2)%d]=domain[j][(axis+2)%d];}
                        group.add(createLine(material,horizontalPoints));

                        for(var j=0;j<2;j++){
                            verticalPoints[j][axis]=domain[w][axis];
                            verticalPoints[j][(axis+2)%d]=domain[0][(axis+2)%d]+i*1.0/segments*edgeLengths[(axis+2)%d];
                            verticalPoints[j][(axis+1)%d]=domain[j][(axis+1)%d];}
                        group.add(createLine(material,verticalPoints));
                    }}}}
        forces.add(group);
    };

    return module;
});
