#include "mex.h"
#include <string.h>
#include "importvol.h"

#define NLHS  1
#define USAGE "usage: volume = loadVol(filename)"

static const char *pchFieldNames[16];

static const char *pchContentType[]       = { "gray", "power", "velocity" };
static mxArray    *mxVolume               = NULL;
static int         iNoOfVolumes           = 1;




static void vPutVolume(stVOLUME *pstVolume,eVOLUMETYPE eVolType)
{
    int      piDims[3];
    int      iNDims = 3;
    mxArray *mxPtr;
    double  *pdData;
    BYTE    *pbyData;
    DWORD    dwVolLen;
    int      iVolNumber = 0;
    
    if (eVolType & eGRAY) {
        iVolNumber = 0;
    } 
    else if (eVolType & eANGIO) {
        iVolNumber = 1;
    }
    else if (eVolType & eCFM) {
        iVolNumber = 2;
    }
    else {
        return;
    }
    if (eVolType & eCARTESIAN) {
        mxSetField(mxVolume,iVolNumber,"type",mxCreateString("cartesian"));
        mxSetField(mxVolume,iVolNumber,"content",mxCreateString(pchContentType[iVolNumber]));
        mxPtr     = mxCreateDoubleMatrix(1,1,mxREAL);
        pdData    = mxGetPr(mxPtr);
        *pdData   = pstVolume->dVoxelSizeMeters;
        mxSetField(mxVolume,iVolNumber,"voxelSize",mxPtr);
        piDims[0] = pstVolume->wSamples;
        piDims[1] = pstVolume->wShots;
        piDims[2] = pstVolume->wImages;
        iNDims    = 3;
        mxPtr     = mxCreateNumericArray(iNDims,piDims,mxUINT8_CLASS,mxREAL);
        pbyData   = (BYTE *)mxGetData(mxPtr);
        dwVolLen  = piDims[0]*piDims[1]*piDims[2];
        memcpy(pbyData,pstVolume->pbyData,pstVolume->dwDataLength * sizeof(BYTE));
        mxSetField(mxVolume,iVolNumber,"data",mxPtr);
    }
    return;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int     buflen = 256;
    int           iNFields = 0;
    char          chFileName[256];
    stVOLUMEFILE *pstVolFile = NULL;
    eVOLUMETYPE   eVolType;
    
    if (nrhs != 1 || nlhs > NLHS) {
        mexErrMsgTxt(USAGE);
    }
    
    
    if ( !mxIsChar(prhs[0]) ) {
        mexErrMsgTxt("Input argument 'filename' must be a string.");
    }
    
    mxGetString(prhs[0],chFileName,buflen-1);
    
    if ((pstVolFile = pstOpenVolumeFile(chFileName)) == NULL) {
        mexErrMsgTxt(pchImportVolErrorMsg());
    }
    eVolType = eVolumeType(pstVolFile);
    if (eVolType & eCFM) {
        iNoOfVolumes = 3;
    }
    else if (eVolType & eANGIO) {
        iNoOfVolumes = 2;
    }
    else {
        iNoOfVolumes = 1;
    }
    if (eVolType & eCARTESIAN) {
        pchFieldNames[0] = "type";
        pchFieldNames[1] = "content";
        pchFieldNames[2] = "voxelSize";
        pchFieldNames[3] = "data";
        iNFields = 4;
    } 
    else {
	mexErrMsgTxt( "unsupported volume format" );
    }
    mxVolume  = mxCreateStructMatrix(iNoOfVolumes,1,iNFields,pchFieldNames);
    plhs[0]   = mxVolume;
    vSetPutVolumeHandler(vPutVolume);
    
    
    if (!iLoadVolumeFile(pstVolFile)) {
        mexErrMsgTxt(pchImportVolErrorMsg());
    }
    vCloseVolumeFile(pstVolFile);
}
