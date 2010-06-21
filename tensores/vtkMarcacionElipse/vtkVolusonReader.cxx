/*=========================================================================
 
        Program:        Visualization Toolkit
        Module:         $RCSfile: vtkVolusonReader.cxx,v $
        Language:       C++
        Date:           $Date: 2001/06/28 19:21:10 $
        Version:        $Version 1.0 $
 
        Autor:          Raul San Jose Estepar
                        ETSI Telecomunicacion, Universidad de Valladolid  Campus Miguel Delibes, Camino del Cementerio, s/n
                       e-mail: rjosest@atenea.tel.uva.es

		Modificado por: Lucilio Cordero Grande
		                ETSI Telecomunicacion, Universidad de Valladolid, Campus Miguel Delibes, Camino del Cementerio, s/n
						e-mail: lcorgra@troya.lpi.tel.uva.es
		Objeto de la modificación: Adaptarlo a una versión más reciente de VTK.
		Fecha modificación: 2007/02/27
 
=========================================================================*/    

#include "vtkVolusonReader.h"
#include "importvol.h"
#include "importvol.cxx"
//#include "vtkScalars.h"
#include "vtkUnsignedCharArray.h"
#include "vtkStructuredPoints.h"
#include "vtkPointData.h"
static const char *pchFieldNames[16];
 static const char *pchContentType[]       = { "gray", "power", "velocity" };

static vtkStructuredPoints *output = NULL;
vtkUnsignedCharArray *outscalars= NULL;
static int         iNoOfVolumes           = 1; 

vtkVolusonReader::vtkVolusonReader()
{
this->FileName = NULL;
}

vtkVolusonReader::~vtkVolusonReader()
{
 if (this->FileName)
    {
    delete [] this->FileName;
    }                    
}

unsigned long int vtkVolusonReader::GetMTime()
{
  unsigned long dtime = this->vtkSource::GetMTime();
  unsigned long rtime = this->vtkSource::GetMTime();
  return (dtime > rtime ? dtime : rtime);
}


static void vPutVolume(stVOLUME *pstVolume,eVOLUMETYPE eVolType)
{
    int      piDims[3];
    int      iNDims = 3;
    double  *pdData;
    unsigned char    *pbyData;
    int    dwVolLen;
    int      iVolNumber = 0;
    double spacing;
	
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
        //mxSetField(mxVolume,iVolNumber,"type",mxCreateString("cartesian"));
        //mxSetField(mxVolume,iVolNumber,"content",mxCreateString(pchContentType[iVolNumber]));
        //mxPtr     = mxCreateDoubleMatrix(1,1,mxREAL);
        //pdData    = mxGetPr(mxPtr);
        //*pdData   = pstVolume->dVoxelSizeMeters;
        //mxSetField(mxVolume,iVolNumber,"voxelSize",mxPtr);
        //Set spacing
		spacing=pstVolume->dVoxelSizeMeters;	
		output->SetSpacing(spacing*1000.0,spacing*1000.0,spacing*1000.0);
		
		//Set dimensions
		piDims[0] = pstVolume->wSamples;
        piDims[1] = pstVolume->wShots;
        piDims[2] = pstVolume->wImages;
        output->SetDimensions(piDims);	
		cout<<"Dimensiones vol salida: "<<piDims[0]<<" "<<piDims[1]<<" "<<piDims[2]<<endl;
		
		//Set scalar data
		//iNDims    = 3;
        //mxPtr     = mxCreateNumericArray(iNDims,piDims,mxUINT8_CLASS,mxREAL);
        //pbyData   = (char *)mxGetData(mxPtr);
        
		//memcpy(pbyData,pstVolume->pbyData,pstVolume->dwDataLength * sizeof(char));
        //mxSetField(mxVolume,iVolNumber,"data",mxPtr);
		
		
		//output->SetScalarType(VTK_CHAR);
		dwVolLen  = piDims[0]*piDims[1]*piDims[2];
		//outscalars=output->GetPointData()->GetScalars();        
		
		//outscalars->SetDataTypeToUnsignedChar();
		//outscalars->SetActiveComponent(iVolNumber);
		outscalars->SetNumberOfValues(pstVolume->dwDataLength);
		
		pbyData=(unsigned char *)outscalars->GetVoidPointer(0);

		memcpy(pbyData,pstVolume->pbyData,pstVolume->dwDataLength * sizeof(char));
    }
    return;
}


 
                                                            

void vtkVolusonReader::Execute()
{

const int     buflen = 256;
int           iNFields = 0;
    
stVOLUMEFILE *pstVolFile = NULL;
eVOLUMETYPE   eVolType;        

output = this->GetOutput();


if ((pstVolFile = pstOpenVolumeFile(this->FileName)) == NULL) {
        vtkErrorMacro("Error Reading Voluson File");
    	exit;
	}        

//cout<<"Vol file opened"<<endl;

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

	outscalars= vtkUnsignedCharArray::New();
	//(output->GetPointData()->GetScalars())->SetNumberOfComponents(iNoOfVolumes);
	//outscalars->SetNumberOfComponents(iNoOfVolumes);

cout<<"Number of Volumes: "<<iNoOfVolumes<<endl;
//cout<<"Assinging components"<<endl;

	vSetPutVolumeHandler(vPutVolume);  	

//cout<<"Before Loading Vol File"<<endl;

if (!iLoadVolumeFile(pstVolFile)) {
        //mexErrMsgTxt(pchImportVolErrorMsg());
    	vtkErrorMacro("Error Loading Voluson File");
	}

//cout<<"After Loading Vol File"<<endl;
    vCloseVolumeFile(pstVolFile);
    outscalars->SetName("Intensidades");                                
	output->GetPointData()->SetScalars(outscalars);
	//output->SetActiveScalars("Intensidades");
	//output->Update();
	//output->UpdateData();
}

void vtkVolusonReader::PrintSelf(ostream& os, vtkIndent indent)
{

 
}                                  



	                       
