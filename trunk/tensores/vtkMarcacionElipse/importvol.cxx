#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <windows.h>
//#include <mex.h>
#include "importvol.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPoints.h"

#ifdef DEBUG
#define ERROR_MSG(text) sprintf(pchErrorMsg,"%s(%d): %s",__FILE__,__LINE__,text)
#else
//#define ERROR_MSG(text) sprintf(pchErrorMsg,"%s",text)
#define ERROR_MSG(text) printf("%s\n",text)
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
	int i;
    unsigned int              dwFileLength;
    char               pchFileID[16];
    unsigned int              dwLength, dwRestInBuffer;
    stVOLUMEFILE      *pstVolumeFile = NULL;
    unsigned char              *pbyBuffer;
    FILE              *pFile;
    unsigned char              *pbyPtr;
    unsigned short              *pwPtr;
    unsigned int             *pdwPtr;
    unsigned short               wGroup;
    unsigned short             wTag;
    unsigned int              dwFileType;

	char auxChar;
	char *auxCharPtr;
	unsigned char auxUChar;
	unsigned char *auxUCharPtr;
	short auxShort;
	short *auxShortPtr;
	unsigned short auxUShort;
	unsigned short *auxUShortPtr;
	int auxInt;
	int *auxIntPtr;
	unsigned int auxUInt;
	unsigned int *auxUIntPtr;
	float auxFloat;
	float *auxFloatPtr;
	double auxDouble;
	double *auxDoublePtr;
	float auxVector[3];
	float auxMatriz[9];


	vtkPoints *Puntos=vtkPoints::New();
	vtkPolyData *Datos=vtkPolyData::New();
	vtkPolyDataWriter *Escritor=vtkPolyDataWriter::New();
    

	
	if ((pFile = fopen(pchFileName,"rb")) == NULL) {
        ERROR_MSG("cannot open file");
        goto error;
    }
    

	
	fseek(pFile,0L,SEEK_END);
    

	
	dwFileLength = (unsigned int)ftell(pFile);

    rewind(pFile);
	
	printf("Data area size: %d\n",dwFileLength);;

	//printf("Before Reading heading\n");
	
    if (fread(pchFileID,sizeof(char),16,pFile) != 16) {
        fclose(pFile);
        ERROR_MSG("read error");
        goto error;
    }
	
	//printf("After Reading heading\n");

	
    /*if (strncmp(pchFileID,FILE_ID,16)) {
	ERROR_MSG("incorrect file ID");
        goto error;
    }*/
    dwFileLength -= 16;
    if ((pbyBuffer = (unsigned char *)malloc(dwFileLength * sizeof(unsigned char))) == NULL) {
        ERROR_MSG("memory allocation error");
        goto error;
    }
    /* read file */
    //printf("Before reading File\n");
	if (fread(pbyBuffer,sizeof(unsigned char),dwFileLength,pFile) != dwFileLength) {
        fclose(pFile);
        ERROR_MSG("read error");
        goto error;
    }
	
	//printf("After reading File\n");
    
	fclose(pFile);
    
	//printf("After closing file\n");
	
	
	if ((pstVolumeFile = (stVOLUMEFILE *) calloc(1,sizeof(struct stVOLUMEFILE))) == NULL) {
        ERROR_MSG("memory allocation error");
        goto error;
    }
	
	//printf("Filling the structure\n");
	
    pstVolumeFile->dwLength  = dwFileLength;
    pstVolumeFile->pbyBuffer = pbyBuffer;
    pstVolumeFile->eType     = eNO_OF_VOLUMETYPES;
    /* parse file to get filetype */
    dwRestInBuffer = pstVolumeFile->dwLength;
    pbyPtr         = pstVolumeFile->pbyBuffer;
    
	//printf("Space left in buffer: %d\n",dwRestInBuffer);
	//printf("Parsing volume type\n");
	
	while (1) {

        if (dwRestInBuffer < 8) {
            ERROR_MSG("unexpected EOF: 1");
            goto error;
        }
        pwPtr    = (unsigned short *)pbyPtr;
	wGroup   = *pwPtr++;
	wTag     = *pwPtr++;
        pdwPtr   = (unsigned int *)pwPtr;
	dwLength = *pdwPtr++;
        pbyPtr   = (unsigned char *)pdwPtr;
        dwRestInBuffer -= 8;
       Puntos->InsertNextPoint((float)dwLength,(float)wGroup,(float)wTag);
	   if (dwLength==1)
	   {
		   auxUCharPtr=(unsigned char *)pbyPtr;
		   auxUChar=*auxUCharPtr;
		   auxCharPtr=(char *)pbyPtr;
		   auxChar=*auxCharPtr;
		   Puntos->InsertNextPoint(auxUChar,auxChar,222);
		   Puntos->InsertNextPoint(222,222,222);
	   }
	   else if (dwLength==2)
	   {
		   auxUShortPtr=(unsigned short *)pbyPtr;
		   auxUShort=*auxUShortPtr;
		   auxShortPtr=(short *)pbyPtr;
		   auxShort=*auxShortPtr;
		   Puntos->InsertNextPoint(auxUShort,auxShort,333);
		   Puntos->InsertNextPoint(333,333,333);
	   }
	   else if (dwLength==4)
	   {
		   auxUIntPtr=(unsigned int *)pbyPtr;
		   auxUInt=*auxUIntPtr;
		   auxIntPtr=(int *)pbyPtr;
		   auxInt=*auxIntPtr;
		   auxFloatPtr=(float *)pbyPtr;
		   auxFloat=*auxFloatPtr;
		   Puntos->InsertNextPoint(auxUInt,auxInt,auxFloat);
		   Puntos->InsertNextPoint(444,444,444);
	   }
	   else if (dwLength==8)
	   {
		   auxDoublePtr=(double *)pbyPtr;
		   auxDouble=*auxDoublePtr;
		   Puntos->InsertNextPoint(auxDouble,555,555);
		   Puntos->InsertNextPoint(555,555,555);
	   }
	   else if (dwLength==12)
	   {
		   for (i=0;i<3;i++)
		   {
			auxFloatPtr=(float *)(pbyPtr+i*4);
		    auxVector[i]=*auxFloatPtr;
		   }
		   Puntos->InsertNextPoint(auxVector[0],auxVector[1],auxVector[2]);
		   Puntos->InsertNextPoint(666,666,666);
	   }
	   else if (dwLength==36)
	   {
		   for (i=0;i<9;i++)
		   {
			auxFloatPtr=(float *)(pbyPtr+i*4);
		    auxMatriz[i]=*auxFloatPtr;
		   }
		   Puntos->InsertNextPoint(auxMatriz[0],auxMatriz[1],auxMatriz[2]);
		   Puntos->InsertNextPoint(auxMatriz[3],auxMatriz[4],auxMatriz[5]);
	   }
	   else if (dwLength<1000)
	   {
		   Puntos->InsertNextPoint(777,777,777);
		   Puntos->InsertNextPoint(777,777,777);
	   }
	   else
	   {
		   Puntos->InsertNextPoint(999,999,999);
		   Puntos->InsertNextPoint(999,999,999);
	   }



		//printf("Rest in Buffer: %d\n",dwRestInBuffer);
		//printf("Length: %d\n",dwLength);
		
		if (dwRestInBuffer < dwLength) {
            ERROR_MSG("unexpected EOF: 2");
            goto error;
        }
        /* interpret goups and tags */
	if (wGroup == 0xffff && wTag == 0xffff) {
		Datos->SetPoints(Puntos);
		Datos->Update();
		Escritor->SetInput(Datos);
		Escritor->SetFileName("Kretz.vtk");
		Escritor->Write();
	    break;              /* EOF */
        }
	
        /* Check File Type */
        if (wGroup == 0x0001 && wTag == 0x0100) { /* file type */
            pdwPtr         = (unsigned int *)pbyPtr; 
            dwFileType     = *pdwPtr;
            if (dwFileType & eCARTESIAN) {
                pstVolumeFile->eType = (enum eVOLUMETYPE) dwFileType;
            }
			else if (dwFileType & eCURVILINEO)
				pstVolumeFile->eType = (enum eVOLUMETYPE) dwFileType;
            else {
				Puntos->InsertNextPoint((float)wGroup,(float)wTag,(float)dwLength);
				Datos->SetPoints(Puntos);
				Datos->Update();
				Escritor->SetInput(Datos);
				Escritor->SetFileName("Kretz.vtk");
				Escritor->Write();
                ERROR_MSG("unsupported filetype");
                goto error;
            }
        }
        else if (wGroup == 0x0002 && (wTag == 0x1000 || wTag == 0x1100 || wTag == 0x1200)) {
			Puntos->InsertNextPoint((float)wGroup,(float)wTag,(float)dwLength);
			Datos->SetPoints(Puntos);
			Datos->Update();
			Escritor->SetInput(Datos);
			Escritor->SetFileName("Kretz.vtk");
			Escritor->Write();
            ERROR_MSG("reading of compressed files not supported");
            goto error;
        }
        pbyPtr         += dwLength;
        dwRestInBuffer -= dwLength;
    }
    if (pstVolumeFile->eType == eNO_OF_VOLUMETYPES) {
		Puntos->InsertNextPoint((float)wGroup,(float)wTag,(float)dwLength);
		Datos->SetPoints(Puntos);
		Datos->Update();
		Escritor->SetInput(Datos);
		Escritor->SetFileName("Kretz.vtk");
		Escritor->Write();
        ERROR_MSG("unsupported filetype");
        goto error;
    }
    
	//printf("Before leaving the reading function\n");
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
    unsigned int              dwRestInBuffer;
    unsigned char              *pbyPtr;
    unsigned short              *pwPtr;
    unsigned int             *pdwPtr;
    unsigned short               wGroup;
    unsigned short               wTag;
    unsigned int              dwTagLength;
    unsigned short               wValue;
    double             dValue;
    int                iVolNo;
    struct stVOLUME   *pstCurrentVolume = NULL;

    if (pstVolumeFile == NULL) {
        ERROR_MSG("file not opened");
        goto error;
    }
    if (pstVolumeFile->eType == eNO_OF_VOLUMETYPES) {
        ERROR_MSG("unsupported filetype");
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
        pwPtr       = (unsigned short *)pbyPtr;
	wGroup      = *pwPtr++;
	wTag        = *pwPtr++;
        pdwPtr      = (unsigned int *)pwPtr;
	dwTagLength = *pdwPtr++;
        pbyPtr      = (unsigned char *)pdwPtr;
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
                wValue = *(unsigned short *)pbyPtr;
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
                pstCurrentVolume->pbyData      = (unsigned char *)pbyPtr;
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
        if ((unsigned int)pstVolume->wSamples * pstVolume->wShots * pstVolume->wImages == pstVolume->dwDataLength) {
            return 1;
        }
    }
    return 0;
}


