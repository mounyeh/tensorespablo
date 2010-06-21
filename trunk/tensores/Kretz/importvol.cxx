#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "importvol.h"

namespace itk
{


#ifdef DEBUG
#define ERROR_MSG(text) sprintf(pchErrorMsg,"%s(%d): %s",__FILE__,__LINE__,text)
#else
#define ERROR_MSG(text) sprintf(pchErrorMsg,"%s",text)
#endif


#define VOLUME_NO(wTag) ((wTag & 0x0F00) >> 8)

#define FILE_ID "KRETZFILE 1.0   "

static char pchErrorMsg[256] = { '\0' };
static int  iIsVolume(struct stVOLUME *pstVolume, enum eVOLUMETYPE eType);

static void (*pvPutVolume)       (struct stVOLUME *,enum eVOLUMETYPE)                                     = NULL;


const char *pchImportVolErrorMsg(void)
{
    return pchErrorMsg;
}

void vSetPutVolumeHandler(void (*vPutVolume)(struct stVOLUME *, enum eVOLUMETYPE))
{
    pvPutVolume        = vPutVolume;
    return;
}


stVOLUMEFILE *pstOpenVolumeFile(const char *pchFileName)
{
    DWORD              dwFileLength;
    char               pchFileID[16];
    DWORD              dwLength, dwRestInBuffer;
    stVOLUMEFILE      *pstVolumeFile = NULL;
    BYTE              *pbyBuffer;
    FILE              *pFile;
    BYTE              *pbyPtr;
    WORD              *pwPtr;
    DWORD             *pdwPtr;
    WORD               wGroup;
    WORD               wTag;
    DWORD              dwFileType;
	char               errorSt[200];
    
    if ((pFile = fopen(pchFileName,"rb")) == NULL) {
        ERROR_MSG("cannot open file");
        goto error;
    }
    fseek(pFile,0L,SEEK_END);
    dwFileLength = (DWORD)ftell(pFile);
    rewind(pFile);
    if (fread(pchFileID,sizeof(char),16,pFile) != 16) {
        fclose(pFile);
        ERROR_MSG("read error");
        goto error;
    }
    if (strncmp(pchFileID,FILE_ID,16)) {
	ERROR_MSG("incorrect file ID");
        goto error;
    }
    dwFileLength -= 16;
    if ((pbyBuffer = (BYTE *)malloc(dwFileLength * sizeof(BYTE))) == NULL) {
        ERROR_MSG("memory allocation error");
        goto error;
    }
    /* read file */
    if (fread(pbyBuffer,sizeof(BYTE),dwFileLength,pFile) != dwFileLength) {
        fclose(pFile);
        ERROR_MSG("read error");
        goto error;
    }
    fclose(pFile);
    if ((pstVolumeFile = (itk::stVOLUMEFILE*)calloc(1,sizeof(struct stVOLUMEFILE))) == NULL) {
        ERROR_MSG("memory allocation error");
        goto error;
    }
    pstVolumeFile->dwLength  = dwFileLength;
    pstVolumeFile->pbyBuffer = pbyBuffer;
    pstVolumeFile->eType     = eNO_OF_VOLUMETYPES;
    /* parse file to get filetype */
    dwRestInBuffer = pstVolumeFile->dwLength;
    pbyPtr         = pstVolumeFile->pbyBuffer;
    while (1) {
        if (dwRestInBuffer < 8) {
            ERROR_MSG("unexpected EOF");
            goto error;
        }
        pwPtr    = (WORD *)pbyPtr;
	wGroup   = *pwPtr++;
	wTag     = *pwPtr++;
        pdwPtr   = (DWORD *)pwPtr;
	dwLength = *pdwPtr++;
        pbyPtr   = (BYTE *)pdwPtr;
        dwRestInBuffer -= 8;
        if (dwRestInBuffer < dwLength) {
            ERROR_MSG("unexpected EOF");
            goto error;
        }
        /* interpret goups and tags */
	if (wGroup == 0xffff && wTag == 0xffff) {
	    break;              /* EOF */
        }
        /* Check File Type */
        if (wGroup == 0x0001 && wTag == 0x0100) { /* file type */
            pdwPtr         = (DWORD *)pbyPtr; 
            dwFileType     = *pdwPtr;
            if (dwFileType & eCARTESIAN) {
                pstVolumeFile->eType = (itk::eVOLUMETYPE)dwFileType;
            }
            else {
				sprintf( errorSt, "unsupported filetype (line 118): dwFileType = %i; eCARTESIAN = %i; dwFileType & eCARTESIAN = %i", dwFileType, eCARTESIAN, dwFileType & eCARTESIAN );
                ERROR_MSG( errorSt );
                goto error;
            }
        }
        else if (wGroup == 0x0002 && (wTag == 0x1000 || wTag == 0x1100 || wTag == 0x1200)) {
            ERROR_MSG("reading of compressed files not supported");
            goto error;
        }
        pbyPtr         += dwLength;
        dwRestInBuffer -= dwLength;
    }
    if (pstVolumeFile->eType == eNO_OF_VOLUMETYPES) {
        ERROR_MSG("unsupported filetype (line 130)");
        goto error;
    }
    return pstVolumeFile;
error:
    vCloseVolumeFile(pstVolumeFile);
    return NULL;
}

void vCloseVolumeFile(stVOLUMEFILE *pstVolumeFile)
{
    if (pstVolumeFile) {
        if (pstVolumeFile->pbyBuffer) {
            free(pstVolumeFile->pbyBuffer);
        }
        free(pstVolumeFile);
    }
    return;
}

int
iLoadVolumeFile(stVOLUMEFILE *pstVolumeFile)
{
    DWORD              dwRestInBuffer;
    BYTE              *pbyPtr;
    WORD              *pwPtr;
    DWORD             *pdwPtr;
    WORD               wGroup;
    WORD               wTag;
    DWORD              dwTagLength;
    WORD               wValue;
    double             dValue;
    int                iVolNo;
    struct stVOLUME   *pstCurrentVolume = NULL;

    if (pstVolumeFile == NULL) {
        ERROR_MSG("file not opened");
        goto error;
    }
    if (pstVolumeFile->eType == eNO_OF_VOLUMETYPES) {
        ERROR_MSG("unsupported filetype (line 170)");
        goto error;
    }
    dwRestInBuffer = pstVolumeFile->dwLength;
    pbyPtr         = pstVolumeFile->pbyBuffer;
    /* loop over tags */
    while (1) {
        if (dwRestInBuffer < 8) {
            ERROR_MSG("unexpected EOF");
            goto error;
        }
        pwPtr       = (WORD *)pbyPtr;
		wGroup      = *pwPtr++;
		wTag        = *pwPtr++;
        pdwPtr      = (DWORD *)pwPtr;
		dwTagLength = *pdwPtr++;
        pbyPtr      = (BYTE *)pdwPtr;
        dwRestInBuffer -= 8;
        if (dwRestInBuffer < dwTagLength) {
            ERROR_MSG("unexpected EOF");
            goto error;
        }
        /* interpret goups and tags */
		if (wGroup == 0xffff && wTag == 0xffff) {
			break;              /* EOF */
		}
        iVolNo = VOLUME_NO(wTag);
        switch (wGroup) {
        case 0xC000:            /* volume dimensions */
            if (VOLUME_NO(wTag) < 3) {
                pstCurrentVolume = &(pstVolumeFile->pstVolume[VOLUME_NO(wTag)]);
                wValue = *(WORD *)pbyPtr;
                switch (wTag & 0x000F) {
                case 0x001:
                    pstCurrentVolume->wSamples = wValue;
                    break;
                case 0x002:
                    pstCurrentVolume->wShots   = wValue;
                    break;
                case 0x003:
                    pstCurrentVolume->wImages  = wValue;
                    break;
                }
            }
            break;
        case 0xC100:            /* voxelsize in meters */
            if (VOLUME_NO(wTag) < 3) {
                pstCurrentVolume = &(pstVolumeFile->pstVolume[VOLUME_NO(wTag)]);
                dValue = *(double *)pbyPtr;
                pstCurrentVolume->dVoxelSizeMeters = dValue;
            }
            break;
        case 0xD000:
            if (VOLUME_NO(wTag) < 3) {
                pstCurrentVolume               = &(pstVolumeFile->pstVolume[VOLUME_NO(wTag)]);
                pstCurrentVolume->pbyData      = (BYTE *)pbyPtr;
                pstCurrentVolume->dwDataLength = dwTagLength;
            }
            break;
        }
        pbyPtr         += dwTagLength;
        dwRestInBuffer -= dwTagLength;
    }
    if (pvPutVolume) {
        if (pstVolumeFile->eType & eCARTESIAN) {
            if (pstVolumeFile->eType & eCARTESIAN_GRAY && iIsVolume(&(pstVolumeFile->pstVolume[0]),eCARTESIAN)) {
                (*pvPutVolume)(&(pstVolumeFile->pstVolume[0]),eCARTESIAN_GRAY);
            }
            if (pstVolumeFile->eType & eCARTESIAN_CFM && iIsVolume(&(pstVolumeFile->pstVolume[1]),eCARTESIAN)) {
                (*pvPutVolume)(&(pstVolumeFile->pstVolume[1]),eCARTESIAN_CFM);
            }
            if (pstVolumeFile->eType & eCARTESIAN_ANGIO && iIsVolume(&(pstVolumeFile->pstVolume[2]),eCARTESIAN)) {
                (*pvPutVolume)(&(pstVolumeFile->pstVolume[2]),eCARTESIAN_ANGIO);
            }
        }
    }
    return 1;
error:
    return 0;
}


eVOLUMETYPE eVolumeType(stVOLUMEFILE *pstVolumeFile)
{
    if (pstVolumeFile) {
        return pstVolumeFile->eType;
    }
    return eNO_OF_VOLUMETYPES;
}


static int iIsVolume(struct stVOLUME *pstVolume, enum eVOLUMETYPE eType)
{
    if (pstVolume == NULL) {
        return 0;
    }
    if (eType & eCARTESIAN) {
        if ((DWORD)pstVolume->wSamples * pstVolume->wShots * pstVolume->wImages == pstVolume->dwDataLength) {
            return 1;
        }
    }
    return 0;
}


} // end namespace itk


