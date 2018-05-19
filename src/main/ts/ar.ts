interface ARMatrix {
    row: number;
    clm: number;
    m: number[];
}

interface ICP2DCoordT {
    x: number;
    y: number;
}

interface ICP3DCoordT {
    x: number;
    y: number;
    z: number;
}

//
// Convert a camera parameter structure into an OpenGL projection matrix.
//
export function arglCameraFrustumRH(matXc2U: number[][], videoWidth: number, videoHeight: number, focalmin: number, focalmax: number): number[]
{
	let icpara = newArray2d<number>(3, 4);
    let trans = newArray2d<number>(3, 4);
    let p = newArray2d<number>(3, 3);
    let q = newArray2d<number>(4, 4);
	
    let i: number, j: number;
	
    if (arParamDecompMatf(matXc2U, icpara, trans) < 0) {
        return;
    }
	for (i = 0; i < 4; i++) {
        icpara[1][i] = (videoHeight - 1)*(icpara[2][i]) - icpara[1][i];
    }
	
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            p[i][j] = icpara[i][j] / icpara[2][2];
        }
    }
    q[0][0] = (2.0 * p[0][0] / (videoWidth - 1));
    q[0][1] = (2.0 * p[0][1] / (videoWidth - 1));
    q[0][2] = -((2.0 * p[0][2] / (videoWidth - 1))  - 1.0);
    q[0][3] = 0.0;
	
    q[1][0] = 0.0;
    q[1][1] = -(2.0 * p[1][1] / (videoHeight - 1));
    q[1][2] = -((2.0 * p[1][2] / (videoHeight - 1)) - 1.0);
    q[1][3] = 0.0;
	
    q[2][0] = 0.0;
    q[2][1] = 0.0;
    q[2][2] = (focalmax + focalmin)/(focalmin - focalmax);
    q[2][3] = 2.0 * focalmax * focalmin / (focalmin - focalmax);
	
    q[3][0] = 0.0;
    q[3][1] = 0.0;
    q[3][2] = -1.0;
    q[3][3] = 0.0;

    let m_projection = new Array<number>(16);
	
    for (i = 0; i < 4; i++) { // Row.
		// First 3 columns of the current row.
        for (j = 0; j < 3; j++) { // Column.
            m_projection[i + j*4] = q[i][0] * trans[0][j] +
									q[i][1] * trans[1][j] +
									q[i][2] * trans[2][j];
        }
		// Fourth column of the current row.
        m_projection[i + 3*4] = q[i][0] * trans[0][3] +
								q[i][1] * trans[1][3] +
								q[i][2] * trans[2][3] +
								q[i][3];
    }	
    return m_projection;
}


export function arglCameraViewRH(para: number[][], m_modelview: number[], scale = 1)
{
	m_modelview[0 + 0*4] = para[0][0]; // R1C1
	m_modelview[0 + 1*4] = para[0][1]; // R1C2
	m_modelview[0 + 2*4] = para[0][2];
	m_modelview[0 + 3*4] = para[0][3];
	m_modelview[1 + 0*4] = -para[1][0]; // R2
	m_modelview[1 + 1*4] = -para[1][1];
	m_modelview[1 + 2*4] = -para[1][2];
	m_modelview[1 + 3*4] = -para[1][3];
	m_modelview[2 + 0*4] = -para[2][0]; // R3
	m_modelview[2 + 1*4] = -para[2][1];
	m_modelview[2 + 2*4] = -para[2][2];
	m_modelview[2 + 3*4] = -para[2][3];
	m_modelview[3 + 0*4] = 0.0;
	m_modelview[3 + 1*4] = 0.0;
	m_modelview[3 + 2*4] = 0.0;
	m_modelview[3 + 3*4] = 1.0;
	if (scale != 0.0) {
		m_modelview[12] *= scale;
		m_modelview[13] *= scale;
		m_modelview[14] *= scale;
	}
}


function arParamDecompMatf( source: number[][], cpara: number[][], trans: number[][] ): number
{
    let r: number, c: number;
    let Cpara: number[][] = newArray2d<number>(3, 4);
    let rem1: number, rem2: number, rem3: number;
    
    if( source[2][3] >= 0.0 ) {
        for( r = 0; r < 3; r++ ){
            for( c = 0; c < 4; c++ ){
                Cpara[r][c] = source[r][c];
            }
        }
    }
    else {
        for( r = 0; r < 3; r++ ){
            for( c = 0; c < 4; c++ ){
                Cpara[r][c] = -(source[r][c]);
            }
        }
    }
    
    for( r = 0; r < 3; r++ ){
        for( c = 0; c < 4; c++ ){
            cpara[r][c] = 0.0;
        }
    }
    cpara[2][2] = normf( Cpara[2][0], Cpara[2][1], Cpara[2][2] );
    trans[2][0] = Cpara[2][0] / cpara[2][2];
    trans[2][1] = Cpara[2][1] / cpara[2][2];
    trans[2][2] = Cpara[2][2] / cpara[2][2];
    trans[2][3] = Cpara[2][3] / cpara[2][2];
	
    cpara[1][2] = dotf( trans[2][0], trans[2][1], trans[2][2],
                      Cpara[1][0], Cpara[1][1], Cpara[1][2] );
    rem1 = Cpara[1][0] - cpara[1][2] * trans[2][0];
    rem2 = Cpara[1][1] - cpara[1][2] * trans[2][1];
    rem3 = Cpara[1][2] - cpara[1][2] * trans[2][2];
    cpara[1][1] = normf( rem1, rem2, rem3 );
    trans[1][0] = rem1 / cpara[1][1];
    trans[1][1] = rem2 / cpara[1][1];
    trans[1][2] = rem3 / cpara[1][1];
    
    cpara[0][2] = dotf( trans[2][0], trans[2][1], trans[2][2],
                      Cpara[0][0], Cpara[0][1], Cpara[0][2] );
    cpara[0][1] = dotf( trans[1][0], trans[1][1], trans[1][2],
                      Cpara[0][0], Cpara[0][1], Cpara[0][2] );
    rem1 = Cpara[0][0] - cpara[0][1]*trans[1][0] - cpara[0][2]*trans[2][0];
    rem2 = Cpara[0][1] - cpara[0][1]*trans[1][1] - cpara[0][2]*trans[2][1];
    rem3 = Cpara[0][2] - cpara[0][1]*trans[1][2] - cpara[0][2]*trans[2][2];
    cpara[0][0] = normf( rem1, rem2, rem3 );
    trans[0][0] = rem1 / cpara[0][0];
    trans[0][1] = rem2 / cpara[0][0];
    trans[0][2] = rem3 / cpara[0][0];
    
    trans[1][3] = (Cpara[1][3] - cpara[1][2]*trans[2][3]) / cpara[1][1];
    trans[0][3] = (Cpara[0][3] - cpara[0][1]*trans[1][3]
                   - cpara[0][2]*trans[2][3]) / cpara[0][0];
    
    for( r = 0; r < 3; r++ ){
        for( c = 0; c < 3; c++ ){
            cpara[r][c] /= cpara[2][2];
        }
    }
    
    return 0;
}

function normf( a: number, b: number, c: number ): number
{
    return( Math.sqrt( a*a + b*b + c*c ) );
}

function dotf( a1: number, a2: number, a3: number,b1: number, b2: number, b3: number ): number
{
    return( a1 * b1 + a2 * b2 + a3 * b3 );
}


export function arCreateCameraMatrix(xsize: number, ysize: number, FOVy: number): number[][] {

    // looks a lot like a perspective matrix
    let f = ysize/2.0 / Math.tan(FOVy / 2.0);
        
    let matXc2U: number[][] = 
    [
        [
            f, //param->mat[0][0] =   f;
            0, //param->mat[0][1] =   0.0;
            xsize/2, //param->mat[0][2] = xsize/2.0;
            0  //param->mat[0][3] =   0.0;
        ],[
            0, //param->mat[1][0] =   0.0;
            f, //param->mat[1][1] =   f;
            ysize/2, //param->mat[1][2] = ysize/2.0;
            0  //param->mat[1][3] =   0.0;
        ], [
            0, //param->mat[2][0] =   0.0;
            0, //param->mat[2][1] =   0.0;
            1, //param->mat[2][2] =   1.0;
            0  //param->mat[2][3] =   0.0;
        ]
    ];

    return matXc2U;
    
}

export function arGetTransMatSquare(matXc2U: number[][], width: number, screenCoords: ICP2DCoordT[]): number[][] {
    let worldCoords: ICP3DCoordT[] = [
        { 
            x: -width/2, 
            y: width/2, 
            z: 0
        }, 
        {
            x: width/2, 
            y: width/2, 
            z: 0
        },
        {
            x: width/2, 
            y: -width/2, 
            z: 0
        }, 
        {
            x: -width/2, 
            y: -width/2, 
            z: 0   
        }
    ];
    
    let initMatXw2Xc = icpGetInitXw2Xc_from_PlanarData(matXc2U, screenCoords, worldCoords);

    let result = icpPoint( 
        matXc2U, 
        screenCoords, 
        worldCoords, 
        initMatXw2Xc
    );

    return result;
}

function icpPoint(
    patternMatXc2U: number[][], 
    screenCoords: ICP2DCoordT[], 
    worldCoords: ICP3DCoordT[], 
    initMatXw2Xc: number[][], 
    breakLoopErrorThresh= 0.1, 
    breakLoopErrorThresh2 = 4, 
    breakLoopErrorRatioThresh = 0.99,
    maxLoop = 10

): number[][] {
    let U: ICP2DCoordT = {
        x: 0, 
        y: 0
    };
    let J_U_S: number[][];
    let dU: number[], dx: number, dy: number;
    let matXw2U: number[][] = [];
    let matXw2Xc: number[][] = [];
    let dS: number[] = new Array<number>(6);
    let err0: number, err1: number;
    let i: number, j: number;

    let num = Math.min(worldCoords.length, screenCoords.length);

    if( num < 3 ) throw 'not enough points';

    J_U_S = newArray2d<number>(num*2, 6);
    dU = new Array<number>(2 * num);

    // populate matrices
    for( j = 0; j < 3; j++ ) {
        matXw2Xc.push(new Array<number>(4));
        matXw2U.push(new Array<number>(4));
        for( i = 0; i < 4; i++ ) {
            matXw2Xc[j][i] = initMatXw2Xc[j][i];
            //matXw2U[j][i] = initMatXw2Xc[j][i];
        } 
    }

    for( i = 0;; i++ ) {
        // I think this isn't required...
        arUtilMatMul( patternMatXc2U, matXw2Xc, matXw2U );

        err1 = 0.0;
        for( j = 0; j < num; j++ ) {
            if( icpGetU_from_X_by_MatX2U( U, matXw2U, worldCoords[j] ) < 0 ) {
                throw 'bad coordinate or something';
            }
            dx = screenCoords[j].x - U.x;
            dy = screenCoords[j].y - U.y;
            err1 += dx*dx + dy*dy;
            dU[j*2+0] = dx;
            dU[j*2+1] = dy;
        }
        err1 /= num;

        if( err1 < breakLoopErrorThresh ) break;
        if( i > 0 && err1 < breakLoopErrorThresh2 && err1/err0 > breakLoopErrorRatioThresh ) break;
        if( i == maxLoop ) break;
        err0 = err1;

        for( j = 0; j < num; j++ ) {
            if( icpGetJ_U_S( J_U_S, 2*j, patternMatXc2U, matXw2Xc, worldCoords[j] ) < 0 ) {
                throw 'unable to get JUS';
            }
        }
        if( icpGetDeltaS( dS, dU, J_U_S, num*2 ) < 0 ) {
            throw 'unable to get deltaS';
        }

        icpUpdateMat( matXw2Xc, dS );
    }


    //*err = err1;

    return matXw2Xc;
}

function arUtilMatMul( s1: number[][], s2: number[][], d: number[][] ): number
{
    let i: number, j: number;

    for( j = 0; j < 3; j++ ) {
        for( i = 0; i < 4; i++) {
            d[j][i] = s1[j][0] * s2[0][i]
                    + s1[j][1] * s2[1][i]
                    + s1[j][2] * s2[2][i];
        }
        d[j][3] += s1[j][3];
    }

    return 0;
}

function icpGetU_from_X_by_MatX2U( u: ICP2DCoordT, matX2U: number[][], coord3d: ICP3DCoordT ): number
{
    let hx: number, hy: number, h: number;

    hx = matX2U[0][0] * coord3d.x + matX2U[0][1] * coord3d.y
       + matX2U[0][2] * coord3d.z + matX2U[0][3];
    hy = matX2U[1][0] * coord3d.x + matX2U[1][1] * coord3d.y
       + matX2U[1][2] * coord3d.z + matX2U[1][3];
    h  = matX2U[2][0] * coord3d.x + matX2U[2][1] * coord3d.y
       + matX2U[2][2] * coord3d.z + matX2U[2][3];

    if( h == 0.0 ) return -1;

    u.x = hx / h;
    u.y = hy / h;

    return 0;
}

function icpGetJ_U_S( J_U_S: number[][], J_U_S_rowOffset: number, matXc2U: number[][], matXw2Xc: number[][], worldCoord: ICP3DCoordT ): number
{
    let J_Xc_S = newArray2d<number>(3, 6);
    let J_U_Xc = newArray2d<number>(2, 3);
    let Xc: ICP3DCoordT = {
        x: 0, 
        y: 0, 
        z: 0
    };
    let i: number, j: number, k: number;

    if( icpGetJ_Xc_S( J_Xc_S, Xc, matXw2Xc, worldCoord ) < 0 ) {

        return -1;
    }

    if( icpGetJ_U_Xc( J_U_Xc, matXc2U, Xc ) < 0 ) {

        return -1;
    }

    for( j = 0; j < 2; j++ ) {
        for( i = 0; i < 6; i++ ) {
            J_U_S[J_U_S_rowOffset + j][i] = 0.0;
            for( k = 0; k < 3; k++ ) {
                J_U_S[J_U_S_rowOffset + j][i] += J_U_Xc[j][k] * J_Xc_S[k][i];
            }
        }
    }

    return 0;
}

function icpGetJ_Xc_S( J_Xc_S: number[][], cameraCoord: ICP3DCoordT, T0: number[][], worldCoord: ICP3DCoordT ): number
{
    let J_Xc_T = newArray2d<number>(3, 12);
    let J_T_S = newArray2d<number>(12, 6);
    let i: number, j: number, k: number;

    if( worldCoord == null ) {
        console.log('oh no!!');
    }

    cameraCoord.x = T0[0][0]*worldCoord.x + T0[0][1]*worldCoord.y + T0[0][2]*worldCoord.z + T0[0][3];
    cameraCoord.y = T0[1][0]*worldCoord.x + T0[1][1]*worldCoord.y + T0[1][2]*worldCoord.z + T0[1][3];
    cameraCoord.z = T0[2][0]*worldCoord.x + T0[2][1]*worldCoord.y + T0[2][2]*worldCoord.z + T0[2][3];

    J_Xc_T[0][0] = T0[0][0] * worldCoord.x;
    J_Xc_T[0][1] = T0[0][0] * worldCoord.y;
    J_Xc_T[0][2] = T0[0][0] * worldCoord.z;
    J_Xc_T[0][3] = T0[0][1] * worldCoord.x;
    J_Xc_T[0][4] = T0[0][1] * worldCoord.y;
    J_Xc_T[0][5] = T0[0][1] * worldCoord.z;
    J_Xc_T[0][6] = T0[0][2] * worldCoord.x;
    J_Xc_T[0][7] = T0[0][2] * worldCoord.y;
    J_Xc_T[0][8] = T0[0][2] * worldCoord.z;
    J_Xc_T[0][9] = T0[0][0];
    J_Xc_T[0][10] = T0[0][1];
    J_Xc_T[0][11] = T0[0][2];

    J_Xc_T[1][0] = T0[1][0] * worldCoord.x;
    J_Xc_T[1][1] = T0[1][0] * worldCoord.y;
    J_Xc_T[1][2] = T0[1][0] * worldCoord.z;
    J_Xc_T[1][3] = T0[1][1] * worldCoord.x;
    J_Xc_T[1][4] = T0[1][1] * worldCoord.y;
    J_Xc_T[1][5] = T0[1][1] * worldCoord.z;
    J_Xc_T[1][6] = T0[1][2] * worldCoord.x;
    J_Xc_T[1][7] = T0[1][2] * worldCoord.y;
    J_Xc_T[1][8] = T0[1][2] * worldCoord.z;
    J_Xc_T[1][9] = T0[1][0];
    J_Xc_T[1][10] = T0[1][1];
    J_Xc_T[1][11] = T0[1][2];

    J_Xc_T[2][0] = T0[2][0] * worldCoord.x;
    J_Xc_T[2][1] = T0[2][0] * worldCoord.y;
    J_Xc_T[2][2] = T0[2][0] * worldCoord.z;
    J_Xc_T[2][3] = T0[2][1] * worldCoord.x;
    J_Xc_T[2][4] = T0[2][1] * worldCoord.y;
    J_Xc_T[2][5] = T0[2][1] * worldCoord.z;
    J_Xc_T[2][6] = T0[2][2] * worldCoord.x;
    J_Xc_T[2][7] = T0[2][2] * worldCoord.y;
    J_Xc_T[2][8] = T0[2][2] * worldCoord.z;
    J_Xc_T[2][9] = T0[2][0];
    J_Xc_T[2][10] = T0[2][1];
    J_Xc_T[2][11] = T0[2][2];

    if( icpGetJ_T_S( J_T_S ) < 0 ) return -1;

    for( j = 0; j < 3; j++ ) {
        for( i = 0; i < 6; i++ ) {
            J_Xc_S[j][i] = 0.0;
            for( k = 0; k < 12; k++ ) {
                J_Xc_S[j][i] += J_Xc_T[j][k] * J_T_S[k][i];
            }
        }
    }

    return 0;
}

function icpGetJ_T_S( J_T_S: number[][] ): number
{
    J_T_S[0][0] = 0.0;
    J_T_S[0][1] = 0.0;
    J_T_S[0][2] = 0.0;
    J_T_S[0][3] = 0.0;
    J_T_S[0][4] = 0.0;
    J_T_S[0][5] = 0.0;

    J_T_S[1][0] = 0.0;
    J_T_S[1][1] = 0.0;
    J_T_S[1][2] = -1.0;
    J_T_S[1][3] = 0.0;
    J_T_S[1][4] = 0.0;
    J_T_S[1][5] = 0.0;

    J_T_S[2][0] = 0.0;
    J_T_S[2][1] = 1.0;
    J_T_S[2][2] = 0.0;
    J_T_S[2][3] = 0.0;
    J_T_S[2][4] = 0.0;
    J_T_S[2][5] = 0.0;

    J_T_S[3][0] = 0.0;
    J_T_S[3][1] = 0.0;
    J_T_S[3][2] = 1.0;
    J_T_S[3][3] = 0.0;
    J_T_S[3][4] = 0.0;
    J_T_S[3][5] = 0.0;

    J_T_S[4][0] = 0.0;
    J_T_S[4][1] = 0.0;
    J_T_S[4][2] = 0.0;
    J_T_S[4][3] = 0.0;
    J_T_S[4][4] = 0.0;
    J_T_S[4][5] = 0.0;

    J_T_S[5][0] = -1.0;
    J_T_S[5][1] = 0.0;
    J_T_S[5][2] = 0.0;
    J_T_S[5][3] = 0.0;
    J_T_S[5][4] = 0.0;
    J_T_S[5][5] = 0.0;

    J_T_S[6][0] = 0.0;
    J_T_S[6][1] = -1.0;
    J_T_S[6][2] = 0.0;
    J_T_S[6][3] = 0.0;
    J_T_S[6][4] = 0.0;
    J_T_S[6][5] = 0.0;

    J_T_S[7][0] = 1.0;
    J_T_S[7][1] = 0.0;
    J_T_S[7][2] = 0.0;
    J_T_S[7][3] = 0.0;
    J_T_S[7][4] = 0.0;
    J_T_S[7][5] = 0.0;

    J_T_S[8][0] = 0.0;
    J_T_S[8][1] = 0.0;
    J_T_S[8][2] = 0.0;
    J_T_S[8][3] = 0.0;
    J_T_S[8][4] = 0.0;
    J_T_S[8][5] = 0.0;

    J_T_S[9][0] = 0.0;
    J_T_S[9][1] = 0.0;
    J_T_S[9][2] = 0.0;
    J_T_S[9][3] = 1.0;
    J_T_S[9][4] = 0.0;
    J_T_S[9][5] = 0.0;

    J_T_S[10][0] = 0.0;
    J_T_S[10][1] = 0.0;
    J_T_S[10][2] = 0.0;
    J_T_S[10][3] = 0.0;
    J_T_S[10][4] = 1.0;
    J_T_S[10][5] = 0.0;

    J_T_S[11][0] = 0.0;
    J_T_S[11][1] = 0.0;
    J_T_S[11][2] = 0.0;
    J_T_S[11][3] = 0.0;
    J_T_S[11][4] = 0.0;
    J_T_S[11][5] = 1.0;

    return 0;
}

function icpGetJ_U_Xc( J_U_Xc: number[][], matXc2U: number[][], cameraCoord: ICP3DCoordT ): number
{
    let w1: number, w2: number, w3: number, w3_w3: number;

    w1 = matXc2U[0][0] * cameraCoord.x + matXc2U[0][1] * cameraCoord.y + matXc2U[0][2] * cameraCoord.z + matXc2U[0][3];
    w2 = matXc2U[1][0] * cameraCoord.x + matXc2U[1][1] * cameraCoord.y + matXc2U[1][2] * cameraCoord.z + matXc2U[1][3];
    w3 = matXc2U[2][0] * cameraCoord.x + matXc2U[2][1] * cameraCoord.y + matXc2U[2][2] * cameraCoord.z + matXc2U[2][3];

    if( w3 == 0.0 ) return -1;

    w3_w3 = w3 * w3;
    J_U_Xc[0][0] = (matXc2U[0][0] * w3 - matXc2U[2][0] * w1) / w3_w3;
    J_U_Xc[0][1] = (matXc2U[0][1] * w3 - matXc2U[2][1] * w1) / w3_w3;
    J_U_Xc[0][2] = (matXc2U[0][2] * w3 - matXc2U[2][2] * w1) / w3_w3;
    J_U_Xc[1][0] = (matXc2U[1][0] * w3 - matXc2U[2][0] * w2) / w3_w3;
    J_U_Xc[1][1] = (matXc2U[1][1] * w3 - matXc2U[2][1] * w2) / w3_w3;
    J_U_Xc[1][2] = (matXc2U[1][2] * w3 - matXc2U[2][2] * w2) / w3_w3;

    return 0;
}

function icpGetDeltaS( S: number[], dU: number[], J_U_S: number[][], n: number ): number
{
    let matS: ARMatrix, matU: ARMatrix, matJ: ARMatrix;
    let matJt: ARMatrix, matJtJ: ARMatrix, matJtU: ARMatrix;

    matS = arMatrixAlloc(6, 1);
    matS.m   = S;

    matU = arMatrixAlloc(n, 1);
    matU.m   = dU;

    matJ = arMatrixAlloc(n, 6);
    //matJ.m   = J_U_S;
    for( let row = 0; row<n; row++ ) {
        for( let col = 0; col<6; col++ ) {
            matJ.m.push(J_U_S[row][col]);
        }
    }

    matJt = arMatrixAllocTrans( matJ );
    matJtJ = arMatrixAllocMul( matJt, matJ );
    matJtU = arMatrixAllocMul( matJt, matU );
    arMatrixSelfInv(matJtJ)
    arMatrixMul( matS, matJtJ, matJtU );

    return 0;
}

function icpUpdateMat( matXw2Xc: number[][], dS: number[] ): number
{
    let q = new Array<number>(7);
    let mat = newArray2d<number>(3, 4);
    let mat2 = newArray2d<number>(3, 4);
    let i: number, j: number;

    if( icpGetQ_from_S(q, dS) < 0 ) return -1;
    if( icpGetMat_from_Q( mat, q ) < 0 ) return -1;

    for( j = 0; j < 3; j++ ) {
        for( i = 0; i < 4; i++ ) {
            mat2[j][i] = matXw2Xc[j][0] * mat[0][i]
                       + matXw2Xc[j][1] * mat[1][i]
                       + matXw2Xc[j][2] * mat[2][i];
        }
        mat2[j][3] += matXw2Xc[j][3];
    }

    for( j = 0; j < 3; j++ ) {
        for( i = 0; i < 4; i++ ) matXw2Xc[j][i] = mat2[j][i];
    }

    return 0;
}

function icpGetQ_from_S( q: number[], s: number[] ): number
{
    let ra: number;

    ra = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
    if( ra == 0.0 ) {
        q[0] = 1.0;
        q[1] = 0.0;
        q[2] = 0.0;
        q[3] = 0.0;
    }
    else {
        ra = Math.sqrt(ra);
        q[0] = s[0] / ra;
        q[1] = s[1] / ra;
        q[2] = s[2] / ra;
        q[3] = ra;
    }
    q[4] = s[3];
    q[5] = s[4];
    q[6] = s[5];

    return 0;
}

function icpGetMat_from_Q( mat: number[][], q: number[] ): number
{
    let cra: number, one_cra: number, sra: number;

    cra = Math.cos(q[3]);
    one_cra = 1 - cra;
    sra = Math.sin(q[3]);

    mat[0][0] = q[0]*q[0]*one_cra + cra;
    mat[0][1] = q[0]*q[1]*one_cra - q[2]*sra;
    mat[0][2] = q[0]*q[2]*one_cra + q[1]*sra;
    mat[0][3] = q[4];
    mat[1][0] = q[1]*q[0]*one_cra + q[2]*sra;
    mat[1][1] = q[1]*q[1]*one_cra + cra;
    mat[1][2] = q[1]*q[2]*one_cra - q[0]*sra;
    mat[1][3] = q[5];
    mat[2][0] = q[2]*q[0]*one_cra - q[1]*sra;
    mat[2][1] = q[2]*q[1]*one_cra + q[0]*sra;
    mat[2][2] = q[2]*q[2]*one_cra + cra;
    mat[2][3] = q[6];

    return 0;
}

function icpGetInitXw2Xc_from_PlanarData(matXc2U: number[][], screenCoords: ICP2DCoordT[], worldCoords: ICP3DCoordT[]): number[][] {
    
    let num = Math.min(worldCoords.length, screenCoords.length);

    let matA: ARMatrix = arMatrixAlloc(num*2, 8, false);
    let matB: ARMatrix = arMatrixAlloc(num*2, 1, false);


    for( let i=0; i<num; i++ ) {
        let worldCoord = worldCoords[i];
        let screenCoord = screenCoords[i];
        matA.m.push(
            worldCoord.x, worldCoord.y, 1, 0, 
            0, 0, -(worldCoord.x * screenCoord.x), -(worldCoord.y * screenCoord.x), 
            0, 0, 0, worldCoord.x, 
            worldCoord.y, 1, -(worldCoord.x * screenCoord.y), -(worldCoord.y * screenCoord.y)
        );
        matB.m.push(
            screenCoord.x, 
            screenCoord.y
        );        
    }

    let matAt = arMatrixAllocTrans(matA);
    let matAtA = arMatrixAllocMul(matAt, matA);
    let matAtB = arMatrixAllocMul(matAt, matB);
    arMatrixSelfInv(matAtA);
    let matC = arMatrixAllocMul(matAtA, matAtB);

    let v: number[][] = [[0, 0, 0],[0, 0, 0],[0, 0, 0]];
    let t: number[] = [0, 0, 0];

    v[0][2] =  matC.m[6];
    v[0][1] = (matC.m[3] - matXc2U[1][2] * v[0][2]) / matXc2U[1][1];
    v[0][0] = (matC.m[0] - matXc2U[0][2] * v[0][2] - matXc2U[0][1] * v[0][1]) / matXc2U[0][0];
    v[1][2] =  matC.m[7];
    v[1][1] = (matC.m[4] - matXc2U[1][2] * v[1][2]) / matXc2U[1][1];
    v[1][0] = (matC.m[1] - matXc2U[0][2] * v[1][2] - matXc2U[0][1] * v[1][1]) / matXc2U[0][0];
    t[2]  =  1.0;
    t[1]  = (matC.m[5] - matXc2U[1][2] * t[2]) / matXc2U[1][1];
    t[0]  = (matC.m[2] - matXc2U[0][2] * t[2] - matXc2U[0][1] * t[1]) / matXc2U[0][0];

    let l1 = Math.sqrt( v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2] );
    let l2 = Math.sqrt( v[1][0]*v[1][0] + v[1][1]*v[1][1] + v[1][2]*v[1][2] );
    v[0][0] /= l1;
    v[0][1] /= l1;
    v[0][2] /= l1;
    v[1][0] /= l2;
    v[1][1] /= l2;
    v[1][2] /= l2;
    t[0] /= (l1+l2)/2;
    t[1] /= (l1+l2)/2;
    t[2] /= (l1+l2)/2;
    if( t[2] < 0.0 ) {
        v[0][0] = -v[0][0];
        v[0][1] = -v[0][1];
        v[0][2] = -v[0][2];
        v[1][0] = -v[1][0];
        v[1][1] = -v[1][1];
        v[1][2] = -v[1][2];
        t[0] = -t[0];
        t[1] = -t[1];
        t[2] = -t[2];
    }

    let check_rotation_result = check_rotation( v );
    if( check_rotation_result != 0 ) {
        console.warn('check rotation returned '+check_rotation_result);
    }

    v[2][0] = v[0][1]*v[1][2] - v[0][2]*v[1][1];
    v[2][1] = v[0][2]*v[1][0] - v[0][0]*v[1][2];
    v[2][2] = v[0][0]*v[1][1] - v[0][1]*v[1][0];
    l1 = Math.sqrt( v[2][0]*v[2][0] + v[2][1]*v[2][1] + v[2][2]*v[2][2] );
    v[2][0] /= l1;
    v[2][1] /= l1;
    v[2][2] /= l1;


    let initMatXw2Xc: number[][] = [];
    for( let row = 0; row<3; row++ ) {
        initMatXw2Xc.push(new Array<number>(4));
    }

    /*
    initMatXw2Xc.set(
        v[0][0], //initMatXw2Xc[0][0] = v[0][0];
        v[1][0], //initMatXw2Xc[0][1] = v[1][0];
        v[2][0], //initMatXw2Xc[0][2] = v[2][0];
        t[0],    //initMatXw2Xc[0][3] = t[0];
        v[0][1], //initMatXw2Xc[1][0] = v[0][1];
        v[1][1], //initMatXw2Xc[1][1] = v[1][1];
        v[2][1], //initMatXw2Xc[1][2] = v[2][1];
        t[1],    //initMatXw2Xc[1][3] = t[1];
        v[0][2], //initMatXw2Xc[2][0] = v[0][2];
        v[1][2], //initMatXw2Xc[2][1] = v[1][2];
        v[2][2], //initMatXw2Xc[2][2] = v[2][2];
        t[2],    //initMatXw2Xc[2][3] = t[2];
        // this line was added because THREE doesn't actually support 4x3 matrices
        0, 0, 0, 1
    );
    */

    initMatXw2Xc[0][0] = v[0][0];
    initMatXw2Xc[1][0] = v[0][1];
    initMatXw2Xc[2][0] = v[0][2];
    initMatXw2Xc[0][1] = v[1][0];
    initMatXw2Xc[1][1] = v[1][1];
    initMatXw2Xc[2][1] = v[1][2];
    initMatXw2Xc[0][2] = v[2][0];
    initMatXw2Xc[1][2] = v[2][1];
    initMatXw2Xc[2][2] = v[2][2];
    initMatXw2Xc[0][3] = t[0];
    initMatXw2Xc[1][3] = t[1];
    initMatXw2Xc[2][3] = t[2];

    return initMatXw2Xc;
}

function arMatrixAlloc(rows: number, cols: number, populateElements?: boolean): ARMatrix {
    let arMatrix: ARMatrix = {
        row: rows, 
        clm: cols, 
        m: []
    };
    if( populateElements ) {
        for( let row=0; row<rows; row++ ) {
            for( let col=0; col<cols; col++ ) {
                let v;
                if( row == col ) {
                    v = 1;
                } else {
                    v = 0;
                }
                arMatrix.m.push(v);
            }
        }
    }
    return arMatrix;
}

function check_rotation( rot: number[][] )
{
    let v1: number[] = [0, 0, 0];
    let v2: number[] = [0, 0, 0];
    let v3: number[] = [0, 0, 0];
    let ca: number, cb: number, k1: number, k2: number, k3: number, k4: number;
    let a: number, b: number, c: number, d: number;
    let p1: number, q1: number, r1: number;
    let p2: number, q2: number, r2: number;
    let p3: number, q3: number, r3: number;
    let p4: number, q4: number, r4: number;
    let w: number;
    let e1: number, e2: number, e3: number, e4: number;
    let rotFlag: number;

    v1[0] = rot[0][0];
    v1[1] = rot[0][1];
    v1[2] = rot[0][2];
    v2[0] = rot[1][0];
    v2[1] = rot[1][1];
    v2[2] = rot[1][2];
    v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
    w = Math.sqrt( v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2] );
    if( w == 0 ) return -1;
    v3[0] /= w;
    v3[1] /= w;
    v3[2] /= w;

    cb = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    if( cb < 0 ) cb *= -1;
    ca = (Math.sqrt(cb+1) + Math.sqrt(1-cb)) * .5;

    if( v3[1]*v1[0] - v1[1]*v3[0] != 0.0 ) {
        rotFlag = 0;
    }
    else {
        if( v3[2]*v1[0] - v1[2]*v3[0] != 0 ) {
            w = v1[1]; v1[1] = v1[2]; v1[2] = w;
            w = v3[1]; v3[1] = v3[2]; v3[2] = w;
            rotFlag = 1;
        }
        else {
            w = v1[0]; v1[0] = v1[2]; v1[2] = w;
            w = v3[0]; v3[0] = v3[2]; v3[2] = w;
            rotFlag = 2;
        }
    }
    if( v3[1]*v1[0] - v1[1]*v3[0] == 0 ) return -1;
    k1 = (v1[1]*v3[2] - v3[1]*v1[2]) / (v3[1]*v1[0] - v1[1]*v3[0]);
    k2 = (v3[1] * ca) / (v3[1]*v1[0] - v1[1]*v3[0]);
    k3 = (v1[0]*v3[2] - v3[0]*v1[2]) / (v3[0]*v1[1] - v1[0]*v3[1]);
    k4 = (v3[0] * ca) / (v3[0]*v1[1] - v1[0]*v3[1]);

    a = k1*k1 + k3*k3 + 1;
    b = k1*k2 + k3*k4;
    c = k2*k2 + k4*k4 - 1;

    d = b*b - a*c;
    if( d < 0 ) return -1;
    r1 = (-b + Math.sqrt(d))/a;
    p1 = k1*r1 + k2;
    q1 = k3*r1 + k4;
    r2 = (-b - Math.sqrt(d))/a;
    p2 = k1*r2 + k2;
    q2 = k3*r2 + k4;
    if ( rotFlag == 1 ) {
        w = q1; q1 = r1; r1 = w;
        w = q2; q2 = r2; r2 = w;
        w = v1[1]; v1[1] = v1[2]; v1[2] = w;
        w = v3[1]; v3[1] = v3[2]; v3[2] = w;
    } else if ( rotFlag == 2 ) {
        w = p1; p1 = r1; r1 = w;
        w = p2; p2 = r2; r2 = w;
        w = v1[0]; v1[0] = v1[2]; v1[2] = w;
        w = v3[0]; v3[0] = v3[2]; v3[2] = w;
    }

    if( v3[1]*v2[0] - v2[1]*v3[0] != 0.0 ) {
        rotFlag = 0;
    }
    else {
        if( v3[2]*v2[0] - v2[2]*v3[0] != 0.0 ) {
            w = v2[1]; v2[1] = v2[2]; v2[2] = w;
            w = v3[1]; v3[1] = v3[2]; v3[2] = w;
            rotFlag = 1;
        }
        else {
            w = v2[0]; v2[0] = v2[2]; v2[2] = w;
            w = v3[0]; v3[0] = v3[2]; v3[2] = w;
            rotFlag = 2;
        }
    }
    if( v3[1]*v2[0] - v2[1]*v3[0] == 0 ) return -1;
    k1 = (v2[1]*v3[2] - v3[1]*v2[2]) / (v3[1]*v2[0] - v2[1]*v3[0]);
    k2 = (v3[1] * ca) / (v3[1]*v2[0] - v2[1]*v3[0]);
    k3 = (v2[0]*v3[2] - v3[0]*v2[2]) / (v3[0]*v2[1] - v2[0]*v3[1]);
    k4 = (v3[0] * ca) / (v3[0]*v2[1] - v2[0]*v3[1]);

    a = k1*k1 + k3*k3 + 1;
    b = k1*k2 + k3*k4;
    c = k2*k2 + k4*k4 - 1;

    d = b*b - a*c;
    if( d < 0 ) return -1;
    r3 = (-b + Math.sqrt(d))/a;
    p3 = k1*r3 + k2;
    q3 = k3*r3 + k4;
    r4 = (-b - Math.sqrt(d))/a;
    p4 = k1*r4 + k2;
    q4 = k3*r4 + k4;
    if ( rotFlag == 1 ) {
        w = q3; q3 = r3; r3 = w;
        w = q4; q4 = r4; r4 = w;
        w = v2[1]; v2[1] = v2[2]; v2[2] = w;
        w = v3[1]; v3[1] = v3[2]; v3[2] = w;
    } else if ( rotFlag == 2 ) {
        w = p3; p3 = r3; r3 = w;
        w = p4; p4 = r4; r4 = w;
        w = v2[0]; v2[0] = v2[2]; v2[2] = w;
        w = v3[0]; v3[0] = v3[2]; v3[2] = w;
    }

    e1 = p1*p3+q1*q3+r1*r3; if( e1 < 0 ) e1 = -e1;
    e2 = p1*p4+q1*q4+r1*r4; if( e2 < 0 ) e2 = -e2;
    e3 = p2*p3+q2*q3+r2*r3; if( e3 < 0 ) e3 = -e3;
    e4 = p2*p4+q2*q4+r2*r4; if( e4 < 0 ) e4 = -e4;
    if( e1 < e2 ) {
        if( e1 < e3 ) {
            if( e1 < e4 ) {
                rot[0][0] = p1;
                rot[0][1] = q1;
                rot[0][2] = r1;
                rot[1][0] = p3;
                rot[1][1] = q3;
                rot[1][2] = r3;
            }
            else {
                rot[0][0] = p2;
                rot[0][1] = q2;
                rot[0][2] = r2;
                rot[1][0] = p4;
                rot[1][1] = q4;
                rot[1][2] = r4;
            }
        }
        else {
            if( e3 < e4 ) {
                rot[0][0] = p2;
                rot[0][1] = q2;
                rot[0][2] = r2;
                rot[1][0] = p3;
                rot[1][1] = q3;
                rot[1][2] = r3;
            }
            else {
                rot[0][0] = p2;
                rot[0][1] = q2;
                rot[0][2] = r2;
                rot[1][0] = p4;
                rot[1][1] = q4;
                rot[1][2] = r4;
            }
        }
    }
    else {
        if( e2 < e3 ) {
            if( e2 < e4 ) {
                rot[0][0] = p1;
                rot[0][1] = q1;
                rot[0][2] = r1;
                rot[1][0] = p4;
                rot[1][1] = q4;
                rot[1][2] = r4;
            }
            else {
                rot[0][0] = p2;
                rot[0][1] = q2;
                rot[0][2] = r2;
                rot[1][0] = p4;
                rot[1][1] = q4;
                rot[1][2] = r4;
            }
        }
        else {
            if( e3 < e4 ) {
                rot[0][0] = p2;
                rot[0][1] = q2;
                rot[0][2] = r2;
                rot[1][0] = p3;
                rot[1][1] = q3;
                rot[1][2] = r3;
            }
            else {
                rot[0][0] = p2;
                rot[0][1] = q2;
                rot[0][2] = r2;
                rot[1][0] = p4;
                rot[1][1] = q4;
                rot[1][2] = r4;
            }
        }
    }

    return 0;
}

function arMatrixAllocTrans(m: ARMatrix): ARMatrix {
    

    let result = arMatrixAlloc(m.clm, m.row, true);
    
	for(let r = 0; r < m.row; r++) {
		for(let c = 0; c < m.clm; c++) {
			arMatrixSet(result, c, r, arMatrixGet(m, r, c));
		}
    }
    return result;
}

function arMatrixAllocMul(a: ARMatrix, b: ARMatrix): ARMatrix {
    let dest = arMatrixAlloc(a.row, b.clm, true);
    arMatrixMul(dest, a, b);
    return dest;
}

function arMatrixMul(dest: ARMatrix, a: ARMatrix, b: ARMatrix) {

    let r, c, i;

	if(a.clm != b.row || dest.row != a.row || dest.clm != b.clm) throw 'matrix cannot be multiplied';

	for(r = 0; r < dest.row; r++) {
		for(c = 0; c < dest.clm; c++) {
			let v = 0;
			for(i = 0; i < a.clm; i++) {
				v += arMatrixGet(a, r, i) * arMatrixGet(b, i, c);
            }
            arMatrixSet(dest, r, c, v);
		}
    }
}

function arMatrixSelfInv(m: ARMatrix): ARMatrix {
    let ap = 0;
    let dimen = m.row;
    let rowa = m.row;
    let wap: number, wcp: number, wbp: number;/* work pointer                 */
    let i: number, j: number, n: number, ip: number, nwork: number;
    let nos: number[] = [];

    let p: number, pbuf: number, work: number;

    switch (dimen) {
            case (0): return(null);                 /* check size */
            case (1): m.m[ap] = 1 / (m.m[ap]);
                      return m;                   /* 1 dimension */
    }

    for (n = 0; n < dimen ; n++)
        nos.push(n);

    for (n = 0; n < dimen ; n++) {
        wcp = ap + n * rowa;
        for (i = n, wap = wcp, p = 0, ip = -1; i < dimen; i++, wap += rowa) {
            if ( p < ( pbuf = Math.abs(m.m[wap])) ) {
                p = pbuf;
                ip = i;
            }
        }
        if (p <= 1.0e-10 || -1 == ip) /* EPS == Threshold value */
            return(null);

        nwork = nos[ip];
        nos[ip] = nos[n];
        nos[n] = nwork;

        for (j = 0, wap = ap + ip * rowa, wbp = wcp; j < dimen ; j++) {
            work = m.m[wap];
            m.m[wap++] = m.m[wbp];
            m.m[wbp++] = work;
        }

        for (j = 1, wap = wcp, work = m.m[wcp]; j < dimen; j++, wap++) {
            m.m[wap] = m.m[(wap + 1)] / work;
        }
        m.m[wap] = 1 / work;

        for (i = 0; i < dimen ; i++) {
            if (i != n) {
                    wap = ap + i * rowa;
                for(j = 1, wbp = wcp, work = m.m[wap]; j < dimen; j++, wap++, wbp++) {
                            m.m[wap] = m.m[(wap + 1)] - work * (m.m[wbp]);
                }
                m.m[wap] = -work * (m.m[wbp]);
            }
        }
    } //end: for (n = 0; n < dimen ; n++)

    for(n = 0; n < dimen ; n++) {
            for(j = n; j < dimen ; j++)
                    if( nos[j] == n) break;
            nos[j] = nos[n];
            for (i = 0, wap = ap + j, wbp = ap + n; i < dimen; i++, wap += rowa, wbp += rowa) {
                    work = m.m[wap];
                    m.m[wap] = m.m[wbp];
                    m.m[wbp] = work;
            }
    }
    return(m);
}

function arMatrixSet(m: ARMatrix, row: number, col: number, value: number) {
    m.m[m.clm * row + col] = value;
}

function arMatrixGet(m: ARMatrix, row: number, col: number): number {
    return m.m[m.clm * row + col];
}

function toMathMatrix(from: ARMatrix): number[][] {
    let m: number[][] = [];
    for( let r=0; r<from.row; r++ ) {
        let mrow = [];
        for( let c=0; c<from.clm; c++ ) {
            mrow.push(arMatrixGet(from, r, c));
        }
        m.push(mrow);
    }
    return m;
}

function fromMathMatrix(from: number[][], to: ARMatrix) {
    for( let r=0; r<to.row; r++ ) {
        for( let c=0; c<to.clm; c++ ) {
            arMatrixSet(to, r, c, from[r][c]);
        }        
    }
    
}

function newArray2d<T>(rows: number, cols: number, f?: (row: number, col: number) => T) {
    let result: T[][] = [];
    for( let rowNumber = 0; rowNumber<rows; rowNumber++ ) {
        let row: T[] = [];
        for( let colNumber=0; colNumber<cols; colNumber++ ) {
            let value: T;
            if( f ) {
                value = f(rowNumber, colNumber);
            }
            row.push(value);
        }
        result.push(row);
    }
    return result;
}