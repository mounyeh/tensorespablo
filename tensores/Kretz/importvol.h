#ifndef __IMPORTVOL_H__
#define __IMPORTVOL_H__
#include <sys/types.h>
#include <stdio.h>

namespace itk
{

#ifdef WIN32
typedef unsigned __int8  BYTE;
typedef unsigned __int16 WORD;
typedef unsigned __int32 DWORD;
#else
typedef u_int8_t  BYTE;
typedef u_int16_t WORD;
typedef u_int32_t DWORD;
#endif

typedef enum eSCANTYPE {
    eTORICAL,
    eLINEAR,
    eTRAPEZOID,
    eLINEAR_STEERED,
    eCONE,
    eNO_OF_SCANTYPES                                                        
} eSCANTYPE;

typedef enum eVOLUMETYPE {
    eCARTESIAN_GRAY    = 0x0001, 
    eCARTESIAN_CFM     = 0x0002, 
    eCARTESIAN_ANGIO   = 0x0004, 
    eCARTESIAN         = 0x000f, 
    eGRAY              = 0x0011,
    eCFM               = 0x0022, 
    eANGIO             = 0x0044, 
    eNO_OF_VOLUMETYPES = 0xff00  
} eVOLUMETYPE;


typedef enum eSLICE3D {
    eSLICE,
    e3D
} eSLICE3D;


typedef struct stVOLUME{
    WORD           wSamples;
    WORD           wShots;
    WORD           wImages;
    double         dVoxelSizeMeters;
    BYTE          *pbyData;
    DWORD          dwDataLength;
} stVOLUME;

typedef struct stVOLUMEFILE{
    DWORD              dwLength;
    BYTE              *pbyBuffer;
    enum eVOLUMETYPE   eType;
    struct stVOLUME    pstVolume[3];
} stVOLUMEFILE;


extern stVOLUMEFILE *pstOpenVolumeFile (const char *pchFileName);
extern void vCloseVolumeFile           (stVOLUMEFILE *pstVolumeFile);
extern int  iLoadVolumeFile            (stVOLUMEFILE *pstVolumeFile);
extern void vSetPutVolumeHandler       (void (*vPutVolume)       (struct stVOLUME *, enum eVOLUMETYPE));
extern eVOLUMETYPE eVolumeType(stVOLUMEFILE *pstVolumeFile);
extern const char *pchImportVolErrorMsg(void);


} // end namespace itk

#endif
