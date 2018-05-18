import * as THREE from 'three';
import * as jsQR from 'jsqr';
import { arCreateCameraMatrix, arGetTransMatSquare } from './ar';

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
    let model = new THREE.Mesh(
        new THREE.CubeGeometry(1, 1, 1), 
        new THREE.MeshBasicMaterial({
            color: 0xFF0000
        })
    );
    let camera: THREE.Camera;

    function onResize() {
        camera = new THREE.PerspectiveCamera(
            45, 
            container.clientWidth/container.clientHeight
        );
        camera.matrixAutoUpdate = false;

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

            // assume d = 1 in world coordinates
            let cx = qr.location.bottomLeftFinderPattern.x + dx/2;
            let cy = qr.location.bottomLeftFinderPattern.y + dy/2;

            try {
                let m = arGetTransMatSquare(arCameraProjectionMatrix, 1, [
                    new THREE.Vector2(qr.location.topLeftCorner.x, qr.location.topLeftCorner.y), 
                    new THREE.Vector2(qr.location.topRightCorner.x, qr.location.topRightCorner.y), 
                    new THREE.Vector2(qr.location.bottomRightCorner.x, qr.location.bottomRightCorner.y), 
                    new THREE.Vector2(qr.location.bottomLeftCorner.x, qr.location.bottomLeftCorner.y),
                    // new THREE.Vector2(-10, -10),
                    // new THREE.Vector2(10, -10), 
                    // new THREE.Vector2(10, 10), 
                    // new THREE.Vector2(-10, 10), 
                ]);
                let tm = from3x4ToMatrix4(m);
                // flip y and z
                /*
                tm = new THREE.Matrix4().set(
                    1, 0, 0, 0, 
                    0, 0, 1, 0, 
                    0, -1, 0, 0, 
                    0, 0, 0, 1,
                ).multiply(tm);
                */
                
                // tm = tm.getInverse(tm);
                let v = new THREE.Vector3(0, 0, 0);
                console.log('v', v.applyMatrix4(tm));    
                let r = new THREE.Euler();
                r.setFromRotationMatrix(tm);
                console.log('r', r);

                arCameraTransform = tm;
            } catch( e ) {
                console.error('uh oh!', e);
            }

        } else {
            arCameraTransform = null;
        }
        videoTexture.needsUpdate = true;
    }

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
                    camera.projectionMatrix.fromArray([1.9102363924347978, 0, 0, 0, 0, -2.5377457054523322, 0, 0, 0.013943280545895442, 0.005830389685211879, 1.0000002000000199, 1, 0, 0, -0.00020000002000000202, 0]);
                        
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

    let update = function() {
        requestAnimationFrame(update);

        if( videoTexture != null ) {
            renderAugmentedReality();
        }

        if( arCameraTransform != null ) {
            //camera.projectionMatrix.copy(from3x4ToMatrix4(arCameraProjectionMatrix));
            //camera.position.set(0, 0, 0).applyMatrix4(arCameraTransform);
            //camera.rotation.setFromRotationMatrix(arCameraTransform);
            //camera.lookAt(new THREE.Vector3(0, 0, 0));

            model.matrix.copy(arCameraTransform);
            model.updateMatrix();
        } else {
            //camera.position.set(5, 5, 5);
            //camera.lookAt(new THREE.Vector3(0, 0, 0));
        }
        camera.updateMatrix();



        renderer.render(scene, camera);

    }
    update();
    
    canvas.onclick = function() {
        if( arCameraProjectionMatrix == null ) {
            arCameraProjectionMatrix = arCreateCameraMatrix(320, 240, Math.PI/4);
        }
        let x = arGetTransMatSquare(arCameraProjectionMatrix, 80, [
            {
                x: 116.85714701674901, 
                y: 73.796876432764705
            }, 
            {
                x: 175.69344358953515, 
                y: 63.579204613264743
            }, 
            {
                x: 184.63623852551325, 
                y: 118.51882412432155
            },
            {
                x: 126.99318088433860, 
                y: 127.34088991574521
            }, 
        ])
        console.log(x);
    }
}