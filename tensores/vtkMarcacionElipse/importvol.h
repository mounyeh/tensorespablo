#ifndef __IMPORTVOL_H__
#define __IMPORTVOL_H__

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
	eCURVILINEO		   = 0x0010,
    eGRAY              = 0x0011,
    eCFM               = 0x0022,
    eANGIO             = 0x0044,
    eNO_OF_VOLUMETYPES = 0xff00
} eVOLUMETYPE;


typedef enum eSLICE3D {
    eSLICE,
    e3D
} eSLICE3D;


typedef struct stVOLUME {
    int           wSamples;
    int           wShots;
    int           wImages;
    double         dVoxelSizeMeters;
    unsigned char          *pbyData;
    long int          dwDataLength;
} stVOLUME;

typedef struct stVOLUMEFILE {
    long int              dwLength;
    unsigned char              *pbyBuffer;
    enum eVOLUMETYPE   eType;
    struct stVOLUME    pstVolume[3];
} stVOLUMEFILE;


extern stVOLUMEFILE *pstOpenVolumeFile (const char *pchFileName);
extern void vCloseVolumeFile           (stVOLUMEFILE *pstVolumeFile);
extern int  iLoadVolumeFile            (stVOLUMEFILE *pstVolumeFile);
extern void vSetPutVolumeHandler       (void (*vPutVolume)       (struct stVOLUME *, enum eVOLUMETYPE));
extern eVOLUMETYPE eVolumeType(stVOLUMEFILE *pstVolumeFile);
extern const char *pchImportVolErrorMsg(void);

#endif
