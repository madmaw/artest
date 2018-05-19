import * as THREE from 'three';
import * as jsQR from 'jsqr';
import * as pako from 'pako';
import { arCreateCameraMatrix, arGetTransMatSquare, arglCameraFrustumRH, arglCameraViewRH } from './ar';

function from3x4ToMatrix4(m: number[][]): THREE.Matrix4 {
    return new THREE.Matrix4().set(
        m[0][0],
        m[0][1],
        m[0][2],
        m[0][3],
        m[1][0],
        m[1][1],
        m[1][2],
        m[1][3],
        m[2][0],
        m[2][1],
        m[2][2],
        m[2][3],
        // fill in missing row
        0,
        0,
        0,
        1,
    );
}

window.onload = function() {

    let container= document.getElementById('container');
    let canvas = <HTMLCanvasElement>document.getElementById('canvas');

    let scene = new THREE.Scene();
    let geometry = new THREE.CubeGeometry(1, 1, 1);
    geometry.faces[0].color.setHex(0xFF0000);
    geometry.faces[1].color.setHex(0xFF0000);
    geometry.faces[2].color.setHex(0x0000FF);
    geometry.faces[3].color.setHex(0x0000FF);
    geometry.faces[4].color.setHex(0x00FF00);
    geometry.faces[5].color.setHex(0x00FF00);
    geometry.faces[6].color.setHex(0xFFFF00);
    geometry.faces[7].color.setHex(0xFFFF00);
    geometry.faces[8].color.setHex(0xFF00FF);
    geometry.faces[9].color.setHex(0xFF00FF);
    geometry.faces[10].color.setHex(0x00FFFF);
    geometry.faces[11].color.setHex(0x00FFFF);

    let model: THREE.Object3D = new THREE.Mesh(
        geometry, 
        new THREE.MeshBasicMaterial({
            color: 0xFFFFFF, 
            vertexColors: THREE.FaceColors,
        })
    );

    function parseModel(input: string): void {
        try {
            let deflated = atob(input);
            let jsonString = pako.inflate(deflated, {
                to: 'string'
            });
            let json = JSON.parse(jsonString);
            let loader = new THREE.ObjectLoader();
            let newModel = loader.parse(json);
            if( model != null ) {
                scene.remove(model);
            }
            model = newModel;
            if( newModel != null ) {
                scene.add(model);
            }
        } catch( e ) {
            console.warn('unable to decompress external model', input, e);
        }

    }

    // attempt to parse the current window object (if any)
    if( window.location.hash != null && window.location.hash.length > 1 ) {
        parseModel(window.location.hash.substr(1));
    }

    let camera: THREE.Camera;

    function onResize() {
        camera = new THREE.PerspectiveCamera(
            45, 
            container.clientWidth/container.clientHeight
        );
        camera.matrixAutoUpdate = false;
        camera.position.set(0, 0, 0);
        camera.lookAt(new THREE.Vector3(0, 0, 0));
        camera.updateMatrix();


        canvas.width = container.clientWidth;
        canvas.height = container.clientHeight;
        canvas.setAttribute('style', `width: ${container.clientWidth}px; height: ${container.clientHeight}px;`);
    }
    window.onresize = onResize;
    onResize();

    let renderer = new THREE.WebGLRenderer({
        canvas: canvas
    });


    scene.add(model);

    scene.background = new THREE.Color(0xFFFFFF);

    let videoMirrored: boolean;
    let videoElement: HTMLVideoElement;
    let videoCanvas: HTMLCanvasElement;
    let videoStream: MediaStream;
    let videoTexture: THREE.Texture;

    let arCameraTransform: THREE.Matrix4;
    let arCameraProjectionMatrix: number[][];

    function startAugmentedReality() {
        navigator.getUserMedia(
            {
                video: {
                    facingMode: {
                        ideal: 'environment'
                    }
                }
            }, 
            (stream: MediaStream) => {
                for( let track of stream.getVideoTracks()) {
                    let facingMode: any = track.getCapabilities().facingMode;
                    if( facingMode != null ) {
                        if( facingMode instanceof Array ) {
                            if( facingMode.length > 0 ) {
                                videoMirrored = facingMode.indexOf('environment') < 0;
                            }
                        } else {
                            videoMirrored = facingMode != 'environment';
                        }    
                    }
                }    
                videoStream = stream;
                videoElement = document.createElement('video');
                videoCanvas = document.createElement('canvas');
                videoElement.srcObject = stream;
                videoElement.onerror = videoElement.onended = stopAugmentedReality;
    
    
                videoElement.addEventListener('loadedmetadata', () => {
                    videoCanvas.width = videoElement.videoWidth;
                    videoCanvas.height = videoElement.videoHeight;
    
                    videoTexture = new THREE.Texture(videoCanvas);
                    videoTexture.wrapS = videoTexture.wrapT = THREE.ClampToEdgeWrapping;
                    videoTexture.minFilter = THREE.LinearFilter;
    
                    arCameraProjectionMatrix = arCreateCameraMatrix(videoElement.videoWidth, videoElement.videoHeight, Math.PI/4);
                    // TODO somehow turn the above into the below?!
                    //camera.projectionMatrix.fromArray([1.9102363924347978, 0, 0, 0, 0, -2.5377457054523322, 0, 0, 0.013943280545895442, 0.005830389685211879, 1.0000002000000199, 1, 0, 0, -0.00020000002000000202, 0]);
                    let m = arglCameraFrustumRH(arCameraProjectionMatrix, videoElement.videoWidth, videoElement.videoHeight, 1, 100);
                    camera.projectionMatrix.fromArray(m);
                                
                    renderAugmentedReality();
    
                    // set the background
                    scene.background = videoTexture;
    
                });
    
            }, 
            (error: MediaStreamError) => {
                // oh well
                console.warn('unable to get camera stream', error);
            }
        );    
    }

    function stopAugmentedReality() {
        if( videoStream != null ) {
            let tracks = videoStream.getTracks();
            for( let track of tracks ) {
                if( track.stop ) {
                    track.stop();
                }
            }
            if( videoStream.stop ) {
                videoStream.stop();
            }
            videoStream = null;
        }
        videoElement = null;
        videoCanvas = null;
        if( videoTexture != null ) {
            videoTexture.dispose();
            videoTexture = null;    
        }
        scene.background = new THREE.Color(0x333333);
    }

    function renderAugmentedReality() {
        let context = videoCanvas.getContext('2d');
        context.save();
        if( videoMirrored ) {
            // un-mirror
            context.scale(-1, 1);
            context.translate(-videoCanvas.width, 0);
        }
        context.drawImage(videoElement, 0, 0, videoCanvas.width, videoCanvas.height);

        context.restore();

        // TODO ensure that the texture is centered and has the correct aspect ratio (currently stretches to fit)
        if( videoElement.videoWidth / videoElement.videoHeight > container.clientWidth / container.clientHeight ) {

        } else {

        }

        
        let imageData = context.getImageData(0, 0, videoCanvas.width, videoCanvas.height);
        let qr = jsQR.default(imageData.data, imageData.width, imageData.height);
        
        if( qr != null ) {


            // TODO calculate d after unrotation
            let dx = qr.location.topRightFinderPattern.x - qr.location.bottomLeftFinderPattern.x;
            let dy = qr.location.topRightFinderPattern.y - qr.location.bottomLeftFinderPattern.y;
            
            // get the actual horizontal/vertical side (not the hypotenuse)
            let d = Math.sqrt((dx * dx + dy * dy)/2);
            //let d = qr.dimension;

            // draw in the detected square
            context.strokeStyle = '#f00';
            context.beginPath();
            context.moveTo(qr.location.topLeftCorner.x, qr.location.topLeftCorner.y);
            context.lineTo(qr.location.topRightCorner.x, qr.location.topRightCorner.y);
            context.lineTo(qr.location.bottomRightCorner.x, qr.location.bottomRightCorner.y);
            context.lineTo(qr.location.bottomLeftCorner.x, qr.location.bottomLeftCorner.y);
            context.closePath();
            context.stroke();

            context.strokeStyle = '#ff0';
            context.beginPath();
            context.arc(qr.location.topLeftCorner.x, qr.location.topLeftCorner.y, d, 0, Math.PI*2);
            context.stroke();

            let cx = qr.location.bottomLeftFinderPattern.x + dx/2;
            let cy = qr.location.bottomLeftFinderPattern.y + dy/2;

            try {
                let m = arGetTransMatSquare(arCameraProjectionMatrix, 1, [
                    new THREE.Vector2(qr.location.topLeftCorner.x, qr.location.topLeftCorner.y), 
                    new THREE.Vector2(qr.location.topRightCorner.x, qr.location.topRightCorner.y), 
                    new THREE.Vector2(qr.location.bottomRightCorner.x, qr.location.bottomRightCorner.y), 
                    new THREE.Vector2(qr.location.bottomLeftCorner.x, qr.location.bottomLeftCorner.y),
                ]);
                //let tm = from3x4ToMatrix4(m);
                let into = new Array<number>(16);
                arglCameraViewRH(m, into);
                let tm = new THREE.Matrix4().fromArray(into);
                arCameraTransform = tm;

                let url = qr.data;
                // is it a URL?
                let parser = document.createElement('a');
                parser.href = url;
                // is there a hidden model?
                let hash = parser.hash;
                if( hash != null && hash.length > 1 ) {
                    // attempt to decode model 
                    parseModel(hash.substr(1));
                }
                

            } catch( e ) {
                console.error('uh oh!', e);
            }

        } else {
            arCameraTransform = null;
        }
        videoTexture.needsUpdate = true;
    }

    let update = function() {
        requestAnimationFrame(update);

        if( videoTexture != null ) {
            renderAugmentedReality();
        }

        if( arCameraTransform != null ) {

            // both these things work? I'm not sure how that's possible (the camera one is more desirable through)

            // model.matrixAutoUpdate = false;
            // model.matrix.copy(arCameraTransform);
            // model.matrixWorldNeedsUpdate = true;

            let m = arglCameraFrustumRH(arCameraProjectionMatrix, videoElement.videoWidth, videoElement.videoHeight, 1, 100);
            camera.projectionMatrix.fromArray(m).multiply(arCameraTransform);
    

        }



        renderer.render(scene, camera);

    }
    update();

    startAugmentedReality();

    /*
    // output some models for compession
    let modelInput = `
    {
        "geometries": [
            {
                "uuid": "1",
                "type": "TorusBufferGeometry",
                "radius": 0.4,
                "tube": 0.2,
                "radialSegments": 9,
                "tubularSegments": 9,
                "arc": 6.283185
            }],
        "materials": [
            {
                "uuid": "2",
                "type": "MeshStandardMaterial",
                "color": 16777215,
                "emissive": 65408,
                "depthFunc": 3
            }],
        "object": {
            "uuid": "3",
            "type": "Mesh",
            "name": "Torus",
            "matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
            "geometry": "1",
            "material": "2"
        }
    }    
    `;
    let modelCompressed = pako.deflate(modelInput, {
        to: 'string'
    });
    let modelBase64 = btoa(modelCompressed);
    console.log(modelBase64);

    parseModel(modelBase64);
    */
}


