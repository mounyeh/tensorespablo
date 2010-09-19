
/* =========================================================================

  Program:   UsimagToolConsole
  Language:  C++
  Date:      5-07-2007
  Version:   1.0

  Copyright (c) 2007 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

========================================================================= */

// My includes
#include "UsimagToolConsole.h"
#include "GenericImageToImageFilter.h"
#include "TensorGUI.h"
// itk includes
#include "tensor/itkTractography.h"
#include "tensor/itkComputeFAFilter.h"
#include <itkMinimumMaximumImageCalculator.h>
#include <itkEventObject.h>
#include <itkSpatialObjectReader.h>
#include <itkLandmarkSpatialObject.h>
#include <itkGradientImageFilter.h>
#include <itkSpatialOrientationAdapter.h>
#include <itkOrientImageFilter.h>
#include "tensor/itkDTITensorToSymTensorImageFilter.h"
#include "strain/itkStrainTensorToSymTensorImageFilter.h"
#include <itkSymmetricSecondRankTensor.h>
#include <itkImageToVectorImageFilter.h>
#include <itkRawImageIO.h>
// vtk includes:
#include "vtkITKUtility.h"
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkContourFilter.h>
#include <vtkStripper.h>
#include <vtkActor2D.h>
#include <vtkRungeKutta4.h>
#include <vtkStreamLine.h>
#include <vtkPointData.h>
#include <vtkAssignAttribute.h>
#include <vtkPLOT3DReader.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkImageGradient.h>
#include <vtkIndent.h>
#include <vtkFieldDataToAttributeDataFilter.h>
#include <vtkCardinalSpline.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>
// std includes
#include <iomanip>
#include <string>
#include <vector>


#include <itkVTKImageIO.h>
#include <itkImageSeriesWriter.h>

#include <strain/itkComputeStrainScalars.h>

UsimagToolConsole::UsimagToolConsole()
{	
	//==========================================================================
	// AN-TONIO
	// Syncronise data vectors with their corresponding viewers:
	m_destino->add("New",0,NULL,NULL,0);
	m_VectorData.addChoice( m_op1 );
	m_VectorData.addChoice( m_op2 );
	m_VectorData.addChoice( m_destino );
	m_VectorData.addBrowser( m_dataBrowser );
	m_VectorModelData.addBrowser( m_modeldataBrowser );
	m_VectorTensorData.addBrowser( m_tensordataBrowser );
	m_VectorDWIData.addBrowser( m_DWIdataBrowser );
	m_VectorSTData.addBrowser( m_strainTensorDataBrowser );
	//==========================================================================
	
	
  //progressSlider->Observe( m_MetaImageReader );
  progressSlider->value( 0 );

  ActiveGroup = m_configPreferences;
  
  for (int i=0; i < 3; i++ ) {
    ImageViewer[i]->SetImage(m_InputImage[i]);
    //ImageViewer[i]->ClickSelectCallBack( ClickSelectCallback, (void *)this );
	ImageViewer[i]->ViewDetails(false);
	ImageViewer[i]->ChangeActiveClass(1);
    ImageViewer[i]->Update();
  }
  ImageColorViewer->SetInputImage(m_RGBImage);
  ImageViewer[0]->SetClickedRadius(2);
  m_NImagesCargadas = 0;
	
  MyTensorConsole->SetActiveInput(m_activeinput);
  MyTensorConsole->SetDataBrowser(m_dataBrowser);
  MyTensorConsole->SetModelDataBrowser(m_modeldataBrowser);
  MyTensorConsole->SetTensorDataBrowser(m_tensordataBrowser);
  MyTensorConsole->SetDWIDataBrowser(m_DWIdataBrowser);
  MyTensorConsole->SetStrainTensorDataBrowser(m_strainTensorDataBrowser);
  MyTensorConsole->SetDestino(m_destino);
  MyTensorConsole->SetOp1(m_op1);
  MyTensorConsole->SetOp2(m_op2);
  MyTensorConsole->SetImageViewers(ImageViewer);
  MyTensorConsole->SetSliders(sliceNumberSlider);
  MyTensorConsole->SetImageViewer3D(ImageViewer3D);
  MyTensorConsole->SetImageViewerStrain3D(ImageViewerStrain3D);
  MyTensorConsole->SetVectorData((void*)&m_VectorData);  
  MyTensorConsole->SetVectorModelData((void*)&m_VectorModelData); 
  MyTensorConsole->SetVectorTensorData((void*)&m_VectorTensorData); 
  MyTensorConsole->SetVectorDWIData((void*)&m_VectorDWIData); 
  MyTensorConsole->SetVectorSTData((void*)&m_VectorSTData); 
	
  MyTensorConsole->SetActiveGroup(ActiveGroup); 
  MyTensorConsole->SetConfigIO(m_configIO); 
  MyTensorConsole->SetImageColorViewer(ImageColorViewer); 
  MyTensorConsole->SetImageOverlay(m_ImageOverlay); 
  MyTensorConsole->SetPanelMessage(panelMessage);
  MyTensorConsole->SetOrientationText(m_OrientationText);
  MyTensorConsole->SetOriginText(m_OriginText);
  MyTensorConsole->SetPanelMeasures(panelMeasures);
  MyTensorConsole->SetMdisplay(Mdisplay);	
  MyTensorConsole->SetFibersText(m_fibersText);
  MyTensorConsole->SetPanelROIsMeasures(panelROIsMeasures);
  MyTensorConsole->SetROIsText(m_ROIsText);	
  MyTensorConsole->SetMinScalarRange(min_scalar_range);
  MyTensorConsole->SetMaxScalarRange(max_scalar_range);	
  MyTensorConsole->Show(); 

  MyTensorConsole->SetActiveStrainImage(m_activeStrainImage);

  MySegmentationConsole->SetActiveInput(m_activeinput);
  MySegmentationConsole->SetDataBrowser(m_dataBrowser);
  MySegmentationConsole->SetModelDataBrowser(m_modeldataBrowser);
  MySegmentationConsole->SetDestino(m_destino);
  MySegmentationConsole->SetOp1(m_op1);
  MySegmentationConsole->SetOp2(m_op2);
  MySegmentationConsole->SetImageViewers(ImageViewer);
  MySegmentationConsole->SetSliders(sliceNumberSlider);
  MySegmentationConsole->SetImageViewer3D(ImageViewer3D);
  MySegmentationConsole->SetVectorData((void*)&m_VectorData);  
  MySegmentationConsole->SetVectorModelData((void*)&m_VectorModelData); 
  MySegmentationConsole->SetActiveGroup(ActiveGroup); 
  MySegmentationConsole->SetConfigIO(m_configIO); 
  MySegmentationConsole->SetThreads(static_cast<unsigned int>(Threads->value())); 
  MySegmentationConsole->SetNImagesCargadas(m_NImagesCargadas); 
  MySegmentationConsole->SetImageColorViewer(ImageColorViewer); 
  MySegmentationConsole->SetProgressSliders(progressSlider); 
  MySegmentationConsole->SetProgressCounter(m_ProgressCounter);
  MySegmentationConsole->SetClassValue(Class_value);
  MySegmentationConsole->SetNumPoints(num_points);
  MySegmentationConsole->Show(); 

  MyFilteringConsole->SetActiveInput(m_activeinput);
  MyFilteringConsole->SetDataBrowser(m_dataBrowser);
  MyFilteringConsole->SetModelDataBrowser(m_modeldataBrowser);
  MyFilteringConsole->SetDestino(m_destino);
  MyFilteringConsole->SetOp1(m_op1);
  MyFilteringConsole->SetOp2(m_op2);
  MyFilteringConsole->SetImageViewers(ImageViewer);
  MyFilteringConsole->SetSliders(sliceNumberSlider);
  MyFilteringConsole->SetImageViewer3D(ImageViewer3D);
  MyFilteringConsole->SetVectorData((void*)&m_VectorData);  
  MyFilteringConsole->SetVectorModelData((void*)&m_VectorModelData); 
  MyFilteringConsole->SetActiveGroup(ActiveGroup); 
  MyFilteringConsole->SetConfigIO(m_configIO); 
  MyFilteringConsole->SetThreads(static_cast<unsigned int>(Threads->value())); 
  MyFilteringConsole->SetNImagesCargadas(m_NImagesCargadas); 
  MyFilteringConsole->SetImageColorViewer(ImageColorViewer); 
  MyFilteringConsole->SetProgressSliders(progressSlider); 
  MyFilteringConsole->SetProgressCounter(m_ProgressCounter);
  MyFilteringConsole->SetClassValue(Class_value);
  MyFilteringConsole->SetNumPoints(num_points);
  MyFilteringConsole->Show(); 
  
  MyBasicOpConsole->SetActiveInput(m_activeinput);
  MyBasicOpConsole->SetDataBrowser(m_dataBrowser);
  MyBasicOpConsole->SetModelDataBrowser(m_modeldataBrowser);
  MyBasicOpConsole->SetDestino(m_destino);
  MyBasicOpConsole->SetOp1(m_op1);
  MyBasicOpConsole->SetOp2(m_op2);
  MyBasicOpConsole->SetImageViewers(ImageViewer);
  MyBasicOpConsole->SetSliders(sliceNumberSlider);
  MyBasicOpConsole->SetImageViewer3D(ImageViewer3D);
  MyBasicOpConsole->SetVectorData((void*)&m_VectorData);  
  MyBasicOpConsole->SetVectorModelData((void*)&m_VectorModelData); 
  MyBasicOpConsole->SetActiveGroup(ActiveGroup); 
  MyBasicOpConsole->SetConfigIO(m_configIO); 
  MyBasicOpConsole->SetThreads(static_cast<unsigned int>(Threads->value())); 
  MyBasicOpConsole->SetImageColorViewer(ImageColorViewer); 
  MyBasicOpConsole->SetProgressSliders(progressSlider); 
  MyBasicOpConsole->SetProgressCounter(m_ProgressCounter);
  MyBasicOpConsole->Show(); 

  m_InFilePageFormat->add("MetaHeader",0,NULL,NULL,0);
  m_InFilePageFormat->add("VOL",0,NULL,NULL,0);
  m_InFilePageFormat->add("Raw",0,NULL,NULL,0);
  m_InFilePageFormat->add("JPEG",0,NULL,NULL,0);
  m_InFilePageFormat->add("PNG",0,NULL,NULL,0);
  m_InFilePageFormat->add("TIFF",0,NULL,NULL,0);
  m_InFilePageFormat->add("DICOM",0,NULL,NULL,0);
  m_InFilePageFormat->add("VTK",0,NULL,NULL,0);
  m_InFilePageFormat->add("Nrrd",0,NULL,NULL,0);
  m_InFilePageFormat->add("Nifti",0,NULL,NULL,0);
  m_InFilePageFormat->add("DICOM Series",0,NULL,NULL,0);
  m_InFilePageFormat->add("Color",0,NULL,NULL,0);
  m_InFilePageFormat->add("DWI",0,NULL,NULL,0);
  
  m_FileFormatPattern[0] = "mha,mhd";
  m_FileFormatPattern[1] = "VOL";
  m_FileFormatPattern[2] = "raw*";
  m_FileFormatPattern[3] = "jpg,jpeg";
  m_FileFormatPattern[4] = "png";
  m_FileFormatPattern[5] = "tif,tiff";
  m_FileFormatPattern[6] = "dcm,ima"; 
  m_FileFormatPattern[7] = "vtk";
  m_FileFormatPattern[8] = "nhdr,nrrd";
  m_FileFormatPattern[9] = "nii,nii.gz";


  m_InFilePageTensorFormat->add("TensorVTK",0,NULL,NULL,0);
  m_InFilePageTensorFormat->add("TensorNrrd",0,NULL,NULL,0);

  m_FileTensorFormatPattern[0] = "vtk";
  m_FileTensorFormatPattern[1] = "nhdr,nrrd";

  m_InByteOrder->add("Little Endian",0,NULL,NULL,0);
  m_InByteOrder->add("Big Endian",0,NULL,NULL,0);
  
  m_InPixelType->add("UCHAR",0,NULL,NULL,0);
  m_InPixelType->add("CHAR",0,NULL,NULL,0);
  m_InPixelType->add("USHORT",0,NULL,NULL,0);
  m_InPixelType->add("SHORT",0,NULL,NULL,0);  
  m_InPixelType->add("FLOAT",0,NULL,NULL,0);
  m_InPixelType->add("DOUBLE",0,NULL,NULL,0);
 
  m_tipoProcesado->add("    Scalar Data",0,NULL,NULL,0);
  m_tipoProcesado->add("    Tensor Data",0,NULL,NULL,0);
  m_tipoProcesado->add("    Model Data",0,NULL,NULL,0);
  m_tipoProcesado->add("    DWI Data",0,NULL,NULL,0);
  m_tipoProcesado->add("    Strain Data",0,NULL,NULL,0);
  m_tipoProcesado->add("    Basic Operations",0,NULL,NULL,0);
  m_tipoProcesado->add("    Filtering",0,NULL,NULL,0);
  m_tipoProcesado->add("    Segmentation",0,NULL,NULL,0);
  m_tipoProcesado->add("    Registration",0,NULL,NULL,0);
  m_tipoProcesado->value(0);

  
  m_Orientation[0]->add("X",0,NULL,NULL,0);
  m_Orientation[0]->add("Y",0,NULL,NULL,0);
  m_Orientation[0]->add("Z",0,NULL,NULL,0);

  m_Orientation[1]->add("X",0,NULL,NULL,0);
  m_Orientation[1]->add("Y",0,NULL,NULL,0);
  m_Orientation[1]->add("Z",0,NULL,NULL,0);
 
  m_Orientation[2]->add("X",0,NULL,NULL,0);
  m_Orientation[2]->add("Y",0,NULL,NULL,0);
  m_Orientation[2]->add("Z",0,NULL,NULL,0);
  
  m_Orientation[3]->add("X",0,NULL,NULL,0);
  m_Orientation[3]->add("Y",0,NULL,NULL,0);
  m_Orientation[3]->add("Z",0,NULL,NULL,0);
  
  m_ImageMode->add("IMG_VAL",0,NULL,NULL,0);
  m_ImageMode->add("IMG_INV",0,NULL,NULL,0);
  m_ImageMode->add("IMG_LOG",0,NULL,NULL,0);
  m_ImageMode->add("IMG_DX",0,NULL,NULL,0);
  m_ImageMode->add("IMG_DY",0,NULL,NULL,0);
  m_ImageMode->add("IMG_DZ",0,NULL,NULL,0);
  m_ImageMode->add("IMG_BLEND",0,NULL,NULL,0);
  m_ImageMode->add("IMG_MIP",0,NULL,NULL,0);
  m_ImageMode->add("IMG_COLOR",0,NULL,NULL,0);
  
  m_colorMode->add("FA",0,NULL,NULL,0);
  m_colorMode->add("ROIs",0,NULL,NULL,0);
  m_colorMode->add("Heat",0,NULL,NULL,0);
  m_colorMode->add("Discrete",0,NULL,NULL,0);
	
  ActiveOptionGroup = m_DataGroup;
	
  Fl_Text_Buffer *buff = new Fl_Text_Buffer();
	buff->text("\n\n\n"
			   "                                                 Saturn ver 1.0\n"
			   "\n"
			   "                                                 Developers: \n" 
			   "       Ruben Cardenes, Emma Munoz Moreno, Antonio Tristan Vega\n"
			   "\n"
			   "\n"
			   "                             Laboratory of Image Processing (LPI)\n"
			   "                                   University of Valladolid (UVA)\n"
			   "                                      ETSI Telecomunications \n"
			   "                                                 Valladolid\n" 
			   "                                                    Spain \n"
			   "\n"
			   "\n"
			   "                               Released only for research purposes\n"
			   "                                             All rights reserved\n");	
  m_about_textdisplay->buffer(buff);
	
  m_destino->value(0);
  num_points->value(0);
	
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** This piece of code is to disable the menus and functionalities that physicians do not require */
#ifdef FORPHYSICIANS
	/*
	 // Alternative code; this will cause the menus to disapear:
	 UsimagToolMainMenu->remove(6);
	 UsimagToolMainMenu->remove(5);
	 UsimagToolMainMenu->remove(4);
	 */
	// Disable the menus corresponding to Basic Ops., Filtering, Segmentation and Registration:
	BasicOpSubMenu->flags = ( BasicOpSubMenu->flags | FL_MENU_INACTIVE );
	FilteringSubMenu->flags = ( FilteringSubMenu->flags | FL_MENU_INACTIVE );
	SegmentationSubMenu->flags = ( SegmentationSubMenu->flags | FL_MENU_INACTIVE );
	RegistrationSubMenu->flags = ( RegistrationSubMenu->flags | FL_MENU_INACTIVE );
	// The help is not prepared yet, so I think it is better to hide this item:
	HelpSubMenu->flags = ( HelpSubMenu->flags | FL_MENU_INACTIVE );
	// All items in "Image" submenu seem to be disabled, so I think it is better to disable the whole submenu:
	ImageSubMenu->flags = ( ImageSubMenu->flags | FL_MENU_INACTIVE );
	// Disable the options in DTI menu which are not available at this moment:
	DTIFilterDWIMenuItem->flags = ( DTIFilterDWIMenuItem->flags | FL_MENU_INACTIVE );
	DTIEstimateTensorMenuItem->flags = ( DTIEstimateTensorMenuItem->flags | FL_MENU_INACTIVE );
	// Disable "Open Raw" menu, since physicians will not use this:
	OpenRawMenuItem->flags = ( OpenRawMenuItem->flags | FL_MENU_INACTIVE );
	// Remove the options which will not be available in the "type of processing" menu:

	m_tipoProcesado->remove( 8 );
	m_tipoProcesado->remove( 7 );
	m_tipoProcesado->remove( 6 );
	m_tipoProcesado->remove( 5 );

/* MODIFICADO AL AÑADIR DATOS DE STRAIN
	m_tipoProcesado->remove( 7 );
	m_tipoProcesado->remove( 6 );
	m_tipoProcesado->remove( 5 );
	m_tipoProcesado->remove( 4 );
*/
	m_tipoProcesado->redraw(); // Force the listbox to update its appearance
#endif
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	
	//==========================================================================
	// Set the number of threads to the default number of threads
	Threads->step( 1 );
	Threads->minimum( 1 );
	Threads->value(   MultiThreader::New()->GetNumberOfThreads()   ); // this makes use of the standard ITK muli-threader
	//==========================================================================

  //====================================================
  /** Set interaction style to trackball by default */
  vtkInteractorStyleSwitch *interactor = 
    vtkInteractorStyleSwitch::SafeDownCast( ImageViewer3D->GetInteractor()->GetInteractorStyle() );
  if ( interactor != NULL )   
     interactor->SetCurrentStyleToTrackballCamera();

  interactor = vtkInteractorStyleSwitch::SafeDownCast( ImageViewerStrain3D->GetInteractor()->GetInteractorStyle() );
  if ( interactor != NULL )   
     interactor->SetCurrentStyleToTrackballCamera();
  //====================================================
}

UsimagToolConsole::~UsimagToolConsole()
{
}


void UsimagToolConsole::Quit()
{
	panel->hide();
}


void UsimagToolConsole::OnFileInputChange(int mode)
{
	// Scalar mode
	if (mode == 0){
		const char *text = m_InFilePageBrowser->value();
		if (text != NULL && strlen(text) > 0){
			m_filename = text;
		}
		
		if (text != NULL && strlen(text) > 0){
			int v = DetermineFileFormatFromFileName(text,0);    
			if (v != -1){
				m_InFilePageFormat->value(v+1);
				m_FileType = m_InFilePageFormat->value();
			}
		}
	}
	
	// Tensor mode
	if (mode == 1){
		const char *text = m_InFilePageTensorBrowser->value();
		if (text != NULL && strlen(text) > 0){
			m_filename = text;
		}
		
		if (text != NULL && strlen(text) > 0){
			int v = DetermineFileFormatFromFileName(text,1);    
			if (v != -1){
				m_InFilePageTensorFormat->value(v+1);
				m_FileType = m_InFilePageTensorFormat->value();
			}
		}
	}
	if (mode == 2){
		const char *text = m_InFilePageDWIBrowser->value();
		if(text != NULL && strlen(text) > 0){
			m_filename = text;
		}
	}
}



int UsimagToolConsole::DetermineFileFormatFromFileName(const char *testFile,int mode)
{
  if (mode == 0) {
    // Iterate over the known file types
    for(int i = 0;i < m_FORMAT_COUNT;i++) {
      // Create a matching pattern
      std::string pattern = "*.{" + m_FileFormatPattern[i] + "}";
    
     // Check if the filename matches the pattern
     if(fl_filename_match(testFile,pattern.c_str()))
       return i;
     }
  }
  // Tensor mode:
  if (mode == 1) {
    // Iterate over the known file types
    for(int i = 0;i < m_TENSOR_FORMAT_COUNT;i++) {
      // Create a matching pattern
      std::string pattern = "*.{" + m_FileTensorFormatPattern[i] + "}";
    
     // Check if the filename matches the pattern
     if(fl_filename_match(testFile,pattern.c_str()))
       return i;
     }
  }

  
  // Failed: return illegal pattern
  return -1;
}

void UsimagToolConsole::OnFileFormatChange(int mode)
{
	//Scalar mode
	if (mode == 0){ 
		m_FileType = m_InFilePageFormat->value();
	}
	// Tensor mode
	if (mode == 1){ 
		m_FileType = m_InFilePageTensorFormat->value();
	}
}

void UsimagToolConsole::OnByteOrderChange()
{
	m_ByteOrder = m_InByteOrder->value();
}


//*** Vero
void UsimagToolConsole::OnTipoProcesadoChange()
{
	switch(m_tipoProcesado->value()){
		case 0:	// Scalar Data
		    ActiveOptionGroup->hide();
			m_DataGroup->show();
			ActiveOptionGroup = m_DataGroup;
			break;
		case 1:	// Tensor Data
		    ActiveOptionGroup->hide();
			m_TensorDataGroup->show();
			ActiveOptionGroup = m_TensorDataGroup;			
			break;
		case 2:	// Model Data
		    ActiveOptionGroup->hide();
			m_ModelDataGroup->show();
			ActiveOptionGroup = m_ModelDataGroup;
			break;
		case 3:	// DWI Data
		    ActiveOptionGroup->hide();
			m_DWIDataGroup->show();
			ActiveOptionGroup = m_DWIDataGroup;
			break;
		case 4:	// Strain Tensor Data
		    ActiveOptionGroup->hide();
			m_StrainTensorDataGroup->show();
			ActiveOptionGroup = m_StrainTensorDataGroup;
			break;

// MODIFICADOS LOS NÚMEROS AL AÑADIR DATOS DE STRAIN:
		case 5: // Basic Operations
			ActiveOptionGroup->hide();
			m_BasicOperations->show();
			ActiveOptionGroup = m_BasicOperations;
			break;
		case 6: // Filtering
		    ActiveOptionGroup->hide();
			m_Filtering->show();
			ActiveOptionGroup = m_Filtering;
			break;
		case 7: // Segmentation
		    ActiveOptionGroup->hide();
			m_Segmentation->show();
			ActiveOptionGroup = m_Segmentation;
			break;
		case 8:	// Registration
		    ActiveOptionGroup->hide();
			m_Registration->show();
			ActiveOptionGroup = m_Registration;
			break;

	}
}
//*** Vero

void UsimagToolConsole::OnPixelTypeChange()
{
  switch(m_InPixelType->value()) { 
  case 0: 
    return;
    break;
  case 1:
    m_PixelType = itk::ImageIOBase::UCHAR;
    break;
  case 2:
    m_PixelType = itk::ImageIOBase::CHAR;
    break;
  case 3:
    m_PixelType = itk::ImageIOBase::USHORT;
    break;
  case 4:
    m_PixelType = itk::ImageIOBase::SHORT;
    break;
  case 5:
    m_PixelType = itk::ImageIOBase::FLOAT;
    break;
  case 6:
    m_PixelType = itk::ImageIOBase::DOUBLE;
    break;
  default:
    break;
  }
}



void UsimagToolConsole::Abrir( int filetype )
{
  switch(m_FileType) { 
  case 0: 
    return;
    break;
  case 1:
    LoadMeta();
    break;
  case 2:
    LoadVOL();
    break;
  case 3:
    LoadRaw();
    break;
  case 4:
    LoadJPEG();
    break;
  case 5:
    LoadPNG();
    break;
  case 6:
    LoadTIFF();
    break;
  case 7:
    LoadDicom();
    break;
  case 8:
    LoadVTK();
    break;
  case 9:
    LoadNrrd();
    break;
  case 10:
    LoadNifti();
    break;
  case 11:
    LoadDicomSeries();
    break;
  case 12:
    LoadColor();
    break;
  case 13:
    LoadGenericDWI();
    break;
  default:
    break;
  }  
}



void UsimagToolConsole::AbrirTensor( int filetype )
{
  switch(m_FileType) { 
  case 0: 
    return;
    break;
  case 1:
    LoadTensorVTK();
    break;
  case 2:
    LoadTensorNrrd();
    break;
  default:
    break;
  }
  
}


void UsimagToolConsole::Save( void )
{
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}
	this->Save( m_filename );
}


/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::Save( const char * filename )
{
	int dataId = m_dataBrowser->value()-1;
	if (dataId<0){return;}
	this->ViewData(m_activeinput,dataId);
	
	FloatWriterType::Pointer imageWriter = FloatWriterType::New();
	imageWriter->SetFileName( filename );
	imageWriter->SetInput( m_VectorData[dataId].image );
	// Attempt to Write
	try{
		imageWriter->Write();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
}

void UsimagToolConsole::SaveOverlay(  )
{
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}
	
	UCharWriterType::Pointer imageWriter = UCharWriterType::New();
	imageWriter->SetFileName( m_filename );
	imageWriter->SetInput( ImageViewer[m_activeinput]->GetInputOverlay() );
	// Attempt to Write
	try{
		imageWriter->Write();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::SaveModel( )
{
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}
	int id = m_modeldataBrowser->value()-1;
	if (id<0){return;}
	
	vtkPolyDataWriter* modelWriter = vtkPolyDataWriter::New();
	modelWriter->SetInput( m_VectorModelData[id].data );
	modelWriter->SetFileName( m_filename );
	
	// Attempt to Write
	try{
		modelWriter->Write();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
}

/** Last modified: 3/06/2009 RUBEN */
void UsimagToolConsole::SaveDWI( )
{
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}
	int id = m_DWIdataBrowser->value()-1;
	if (id<0){return;}
	
	// Obtenemos los datos para la cabecera que escribiremos "manualmente"
	float origin[3];
	origin[0] = (float)m_VectorDWIData[id].image->GetOrigin()[0];
	origin[1] = (float)m_VectorDWIData[id].image->GetOrigin()[1];
	origin[2] = (float)m_VectorDWIData[id].image->GetOrigin()[2];
	float B_Value = m_VectorDWIData[id].image->GetB_Value();
	unsigned int num_grads = m_VectorDWIData[id].image->GetNumImages();
	unsigned int Size[3];
	Size[0] = (unsigned int)m_VectorDWIData[id].image->GetLargestPossibleRegion().GetSize()[0];
	Size[1] = (unsigned int)m_VectorDWIData[id].image->GetLargestPossibleRegion().GetSize()[1];
	Size[2] = (unsigned int)m_VectorDWIData[id].image->GetLargestPossibleRegion().GetSize()[2];
	float spacing[3];
	spacing[0] = (float)m_VectorDWIData[id].image->GetSpacing()[0];
	spacing[1] = (float)m_VectorDWIData[id].image->GetSpacing()[1];
	spacing[2] = (float)m_VectorDWIData[id].image->GetSpacing()[2];
	// Pillamos las direcciones:
	typedef	 itk::Vector<double,3>				        DirectionVectorType;
    std::vector< DirectionVectorType> directions;		
    m_VectorDWIData[id].image->GetDiffusionDirections( directions );
	
	// Mostramos datos de cabecera 
	std::cout << "NumImages (Ngrads+baseline): " << num_grads << std::endl;
	std::cout << "B_Value: " << B_Value << std::endl;
	std::cout << "Region size: " << m_VectorDWIData[id].image->GetLargestPossibleRegion().GetSize() << std::endl;
	std::cout << "Origin: " << m_VectorDWIData[id].image->GetOrigin() << std::endl;
	std::cout << "Spacing: " << m_VectorDWIData[id].image->GetSpacing() << std::endl;
	for (unsigned int i=0;i<num_grads-1;i++) {
	  std::cout << "Directions: " << directions[i][0] << " " << directions[i][1] << " " << directions[i][2] << std::endl;
	}
	
	// Metemos los datos en un vector imagen
	//typedef itk::VectorImage<InputPixelType, 3> vectorImageType;
	//vectorImageType::Pointer vectorImage = vectorImageType::New();
	
	//vectorImage->SetVectorLength(m_VectorDWIData[id].image->GetNumImages());
	//vectorImage->SetRegions(m_VectorDWIData[id].image->GetLargestPossibleRegion());
	//vectorImage->SetOrigin(m_VectorDWIData[id].image->GetOrigin());
	//vectorImage->SetSpacing(m_VectorDWIData[id].image->GetSpacing());
	//vectorImage->Allocate();
	
	//typedef itk::ImageToVectorImageFilter< InputImageType > ImageToVectorfilterType;
	//ImageToVectorfilterType::Pointer filter = ImageToVectorfilterType::New();
	InputImageType::Pointer image_aux[num_grads]; 
	
	
	InputImageType::Pointer totalImage	= InputImageType::New(); 
    InputImageType:SizeType size_total;
	size_total[0] = Size[0];
	size_total[1] = Size[1];
	size_total[2] = Size[2]*num_grads;
	InputImageType::RegionType region;
	region.SetSize(size_total);
	totalImage->SetRegions(region);
	totalImage->SetOrigin(m_VectorDWIData[id].image->GetOrigin());
	totalImage->SetSpacing(m_VectorDWIData[id].image->GetSpacing());
	totalImage->Allocate();
	
	typedef itk::ImageRegionIterator< InputImageType> ImageIteratorType;
	typedef itk::ImageRegionConstIterator< InputImageType> ImageConstIteratorType;
	ImageIteratorType it_v(totalImage, totalImage->GetRequestedRegion());
	it_v.GoToBegin();
	
    for (unsigned int i=0;i<num_grads;i++) {
	  image_aux[i]	= InputImageType::New(); 
	  image_aux[i]->SetRegions(m_VectorDWIData[id].image->GetLargestPossibleRegion());
	  image_aux[i]->SetOrigin(m_VectorDWIData[id].image->GetOrigin());
	  image_aux[i]->SetSpacing(m_VectorDWIData[id].image->GetSpacing());
	  image_aux[i]->Allocate();
	  m_VectorDWIData[id].image->GetImageComponent(i, image_aux[i]);
	  m_VectorDWIData[id].image->Update();
	  //filter->SetNthInput( i , image_aux[i] );
	  ImageConstIteratorType it_aux(image_aux[i], image_aux[i]->GetRequestedRegion());
	  for ( it_aux.GoToBegin(); !it_aux.IsAtEnd(); ++it_aux,++it_v )
			it_v.Set( it_aux.Get() );
	}
    //filter->Update();
	
	
	char m_filename_raw[200];
	char m_filename_nhdr[200];
	sprintf(m_filename_nhdr,"%s.nhdr",m_filename);
	sprintf(m_filename_raw,"%s.raw",m_filename);
	
	// Creamos la cabecera:
	char st[200];
	std::ofstream nhdr_file;
	nhdr_file.open(m_filename_nhdr);
	nhdr_file<< "NRRD0005" << std::endl; 
    nhdr_file<< "type: float" << std::endl; 
	nhdr_file<< "dimension: 4"<< std::endl;
    nhdr_file<< "space: right-anterior-superior" << std::endl;
    nhdr_file<< "sizes: " << Size[0] << " "<< Size[1] << " "<< Size[2] << " " << num_grads << std::endl;
    nhdr_file<< "thicknesses: NaN NaN " << spacing[2] << " NaN "  <<std::endl;
	nhdr_file<< "space directions: ("<< spacing[0] << ",0,0) (0,"<< spacing[1] <<",0) (0,0," << spacing[2] <<") none" << std::endl;
    nhdr_file<< "centerings: cell cell cell none"<<std::endl;
    nhdr_file<< "kinds: space space space list"<<std::endl;
    nhdr_file<< "endian: little" <<std::endl;
    nhdr_file<< "encoding: raw" <<std::endl;
    nhdr_file<< "byte skip: -1" <<std::endl;
	nhdr_file<< "space origin: (" << origin[0]<< "," << origin[1] << ","<< origin[2] << ")" <<std::endl;
    nhdr_file<< "measurement frame: (1,0,0) (0,1,0) (0,0,1)" << std::endl;
    nhdr_file<< "data file: " << m_filename_raw << std::endl;
    nhdr_file<< "modality:=DWMRI" << std::endl;
	nhdr_file<< "DWMRI_b-value:="<< B_Value << std::endl;
    nhdr_file<< "DWMRI_gradient_0000:= 0 0 0" << std::endl;
	for (unsigned int i=0;i<num_grads-1;i++) {
		sprintf(st,"DWMRI_gradient_%04d:= %f %f %f\n",i+1,directions[i][0], directions[i][1], directions[i][2]);
		nhdr_file << st;
    }
	nhdr_file.close();
	
	// Escribimos el fichero raw
	//typedef itk::ImageFileWriter< vectorImageType> VectorWriterType;
	
	//VectorWriterType::Pointer vectorWriter = VectorWriterType::New(); 
	//vectorWriter->SetInput( vectorImage );
	//vectorWriter->SetInput( filter->GetOutput() );
	//vectorWriter->SetImageIO(RawImageIO);
	//vectorWriter->SetFileName( m_filename_raw );
	
	RawReaderType::Pointer RawImageIO = RawReaderType::New();
	RawImageIO->SetByteOrderToLittleEndian();
	
	FloatWriterType::Pointer writer = FloatWriterType::New(); 
	writer->SetInput( totalImage );
	writer->SetFileName( m_filename_raw );
	//RawImageIO->SetNumberOfComponents( 1 );
	//RawImageIO->SetNumberOfDimensions( 3 );
	writer->SetImageIO(RawImageIO);
	
	// Attempt to Write
	try{
		writer->Write();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
}


/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::SaveTensor( )
{
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 ){return;}
	
	TensorWriterType::Pointer tensorWriter = TensorWriterType::New();
	tensorWriter->SetFileName( m_filename );
	tensorWriter->SetInput(m_VectorTensorData[dataId].image);
	
	// Attempt to Write
	try{
		tensorWriter->Write();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::SaveTensorVTK( ) {
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 ){return;}
	
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Move these definitions to UsimagToolBase.h!!! */
	typedef itk::Image<itk::SymmetricSecondRankTensor<float, 3>, 3 >				            SymmetricTensorImageType;
	typedef itk::DTITensorToSymTensorImageFilter<TensorImageType, SymmetricTensorImageType>	    CastFilterType;
	typedef CastFilterType::Pointer													            CastFilterPointerType;
	typedef itk::ImageFileWriter<SymmetricTensorImageType>                                      TensorWriterType;
	

	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	TensorWriterType::Pointer   tensorWriter  =  TensorWriterType::New();
	CastFilterPointerType	    caster        =  CastFilterType::New();
	caster->SetInput( m_VectorTensorData[dataId].image );
	tensorWriter->SetInput( caster->GetOutput() );

	tensorWriter->SetFileName(m_filename);
	tensorWriter->Modified();
	try{
		tensorWriter->Write();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}

}

void UsimagToolConsole::SaveStrainTensor( ) {

	typedef itk::StrainTensorToSymTensorImageFilter<STImageType, VectorImage4DType>			CastFilterType;
	typedef CastFilterType::Pointer								 	CastFilterPointerType;

	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}

	int dataId = m_strainTensorDataBrowser->value()-1;
	if( dataId<0 ){return;}

	STImageType::Pointer strainImage = m_VectorSTData[dataId].image;

	STWriterType::Pointer strainWriter = STWriterType::New();
	CastFilterPointerType caster = CastFilterType::New();
	caster->SetInput( strainImage );

	strainWriter->SetInput( caster->GetOutput() );

	strainWriter->Modified();

	// Se generan los nombres de ficheros, añadiendo número y extensión al indicado
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	std::string format = m_filename;
	if (format.rfind('.') != format.npos)
		format.insert(format.rfind('.'),"%02d");

	else format += "%02d.vtk";
	
	int numFiles = strainImage->GetRequestedRegion().GetSize()[3];
	nameGenerator->SetSeriesFormat(format);
	nameGenerator->SetStartIndex(0);
	nameGenerator->SetEndIndex(numFiles-1);
	nameGenerator->SetIncrementIndex(1);

	strainWriter->SetFileNames(nameGenerator->GetFileNames());

	// En primer lugar se escriben los datos de tensor de esfuerzo
	try{
		strainWriter->Write();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}


	DeformWriterType::Pointer   deformWriter  =  DeformWriterType::New();
	deformWriter->SetInput( m_VectorSTData[dataId].deform_image );

	deformWriter->Modified();

	nameGenerator->SetSeriesFormat(format);
	nameGenerator->SetStartIndex(numFiles);
	nameGenerator->SetEndIndex(2*numFiles-1);
	nameGenerator->SetIncrementIndex(1);

	deformWriter->SetFileNames(nameGenerator->GetFileNames());

	// En segundo lugar se escriben los datos de deformación
	try{
		deformWriter->Write();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}


}

void UsimagToolConsole::Load( void )
{
	m_filename = fl_file_chooser("Image Filename","*.{raw,mhd,mha,VOL,jpeg,jpg,tif,tiff,png,vtk,nhdr,flt,vol,vols,volf,ush,chr,char,nii,nii.gz}","");
	if( !m_filename ){return;}
	m_InFilePageBrowser->value( m_filename );
}

void UsimagToolConsole::LoadTensor( void )
{
	m_filename = fl_file_chooser("Image Filename","*.{raw,mhd,mha,vtk,nhdr,flt}","");
	if( !m_filename ){return;}
	m_InFilePageTensorBrowser->value( m_filename );
}

void UsimagToolConsole::LoadDWI( void )
{
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}
	m_InFilePageDWIBrowser->value( m_filename );
}

void UsimagToolConsole::OnProgress(itk::Object *object, const itk::EventObject &event)
{
	// Get the value of the progress
	itk::ProcessObject *po = reinterpret_cast<ProcessObject *>(object);
	float progress = po->GetProgress();
	float max = 1.0;
	if ( m_InputImage[m_activeinput]->GetLargestPossibleRegion().GetSize()[2] == 1)
		max = 2.0;
	m_ProgressCounter->value(100*progress/max);
	
	// Display the filter's progress
	progressSlider->value(100*progress/max);
	
	// Let the UI refresh
	Fl::check();
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadMeta( void )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&UsimagToolConsole::OnProgress);
	
	if ( m_PixelComponents->value() == 1 ){
		FileReaderType::Pointer  metaImageReader = FileReaderType::New();
		/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
		/** TO DO: Check if this piece of code is really necessary. The vectorial counterpart has not an analogous */
		MetaImageIOType::Pointer metaReader      = MetaImageIOType::New();
		metaImageReader->SetImageIO( metaReader );
		/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
		metaImageReader->SetFileName( m_filename );
		metaImageReader->AddObserver( itk::ProgressEvent(), callback );
		metaImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
		// Attempt to read
		try{
			metaImageReader->Update();
		}
		catch( itk::ExceptionObject & e ){
			fl_alert( e.GetDescription() );
			return;
		}
		this->Distribuir( metaImageReader->GetOutput(), m_activeinput, m_filename );
		metaImageReader = NULL;
	}
	else if( m_PixelComponents->value() == 2){
		FileVectorReaderType::Pointer metaVectorImageReader = FileVectorReaderType::New();
		metaVectorImageReader->SetFileName( m_filename );
		metaVectorImageReader->AddObserver( itk::ProgressEvent(), callback );
		metaVectorImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
		// Attempt to read
		try{
			metaVectorImageReader->Update();
		}
		catch( itk::ExceptionObject & e ){
			fl_alert( e.GetDescription() );
			return;
		}
		/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
		// TO DO: We do absolutely nothing with the output of the reader. Is it necessary?
		/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
		metaVectorImageReader = NULL;
	}
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadNifti( void )
{
	if( !m_filename ){return;}
    CommandPointer callback = CommandType::New();  
    callback->SetCallbackFunction(this,&UsimagToolConsole::OnProgress);
	
	FileReaderType::Pointer   niftiImageReader =  FileReaderType::New();
	NiftiImageIOType::Pointer niftiReader      =  NiftiImageIOType::New();
	
	niftiImageReader->SetFileName(m_filename);
	niftiImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	niftiImageReader->AddObserver( itk::ProgressEvent(), callback );
	niftiImageReader->SetImageIO( niftiReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	
    // Attempt to read
    try{
		niftiImageReader->Update();
    }
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
    }
    this->Distribuir( niftiImageReader->GetOutput(), m_activeinput, m_filename );
	niftiImageReader = NULL;
}


/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadVOL( void )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&UsimagToolConsole::OnProgress);
	
	FileReaderType::Pointer genericReader = FileReaderType::New();
	genericReader->SetFileName( m_filename );
	genericReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	genericReader->AddObserver( itk::ProgressEvent(), callback );
	
	// Attempt to read
	try{
		genericReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	this->Distribuir( genericReader->GetOutput(), m_activeinput, m_filename );
	genericReader = NULL;

}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadColor( void )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	ColorReaderType::Pointer colorReader = ColorReaderType::New();
	colorReader->SetFileName( m_filename );
	colorReader->AddObserver( itk::ProgressEvent(), callback );
	colorReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	
	// Attempt to read
	try {
		colorReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	this->DistribuirColor( colorReader->GetOutput(), m_activeinput, m_filename );
	colorReader = NULL;
}



void UsimagToolConsole::LoadRaw( void )
{
	panelDim->show();
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadRaw( const char * filename ) 
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	FileReaderType::Pointer rawImageReader = FileReaderType::New();
	rawImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	RawReaderType::Pointer  rawReader    = RawReaderType::New();
	rawImageReader->SetImageIO( rawReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	rawImageReader->SetFileName( filename );
	rawReader->AddObserver( itk::ProgressEvent(), callback );
	rawReader->SetFileDimensionality( 3 );
	rawReader->SetDimensions( 0, (unsigned int)m_NumberOfPixelsInX->value() );
	rawReader->SetDimensions( 1, (unsigned int)m_NumberOfPixelsInY->value() );
	rawReader->SetDimensions( 2, (unsigned int)m_NumberOfPixelsInZ->value() );
	rawReader->SetComponentType(m_PixelType);
	
	if (m_ByteOrder == 1)
		rawReader->SetByteOrderToLittleEndian();
	else
		rawReader->SetByteOrderToBigEndian();
	
	// Attempt to read
	try{
		rawImageReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	this->Distribuir( rawImageReader->GetOutput(), m_activeinput, filename );
	rawImageReader = NULL;
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadJPEG(  ) {
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	FileReaderType::Pointer JPEGImageReader = FileReaderType::New();
	JPEGImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	JPEGImageIOType::Pointer  JPEGReader    = JPEGImageIOType::New();
	JPEGImageReader->SetImageIO( JPEGReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	JPEGImageReader->SetFileName( m_filename );
	JPEGImageReader->AddObserver( itk::ProgressEvent(), callback );
	try{
		JPEGImageReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	this->Distribuir( JPEGImageReader->GetOutput(), m_activeinput, m_filename );
	JPEGImageReader = NULL;
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadVTK(  )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	FileReaderType::Pointer VTKImageReader = FileReaderType::New();
	VTKImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	VTKImageReader->SetFileName( m_filename );
	VTKImageReader->AddObserver( itk::ProgressEvent(), callback );
	
	try {
		VTKImageReader->Update();
	}
	catch( itk::ExceptionObject & e ) {
		fl_alert( e.GetDescription() );
		return;
	}
	
	this->Distribuir( VTKImageReader->GetOutput(), m_activeinput, m_filename );
	VTKImageReader = NULL;
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadTensorNrrd(  )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&UsimagToolConsole::OnProgress);
	
	FileTensorReaderType::Pointer tensorNrrdReader = FileTensorReaderType::New();
	tensorNrrdReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	tensorNrrdReader->SetFileName( m_filename );
	tensorNrrdReader->AddObserver( itk::ProgressEvent(), callback );
	
	try {
		tensorNrrdReader->Update();
	}
	catch( itk::ExceptionObject & e ) {
		fl_alert( e.GetDescription() );
		return;
	}
	
	this->DistribuirTensor( tensorNrrdReader->GetOutput(), m_activeinput, m_filename );
	tensorNrrdReader = NULL;
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadTensorVTK(  )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	FileTensorReaderType::Pointer tensorImageReader = FileTensorReaderType::New();
	tensorImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	TensorVTKImageIOType::Pointer VTKReader = TensorVTKImageIOType::New(); 
	tensorImageReader->SetImageIO( VTKReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	
	tensorImageReader->SetFileName( m_filename );
	tensorImageReader->AddObserver( itk::ProgressEvent(), callback );
	
	try{
		tensorImageReader->Update();
	}
	catch( itk::ExceptionObject & e ) {
		fl_alert( e.GetDescription() );
		return;
	}
	
	////////////////////////////////////////////////////////////////////////////////////
	// Si cargamos dos volumenes se ajusta el origen para que tenga el mismo valor en las dos imágenes
	// y así poder compararlos mejor visualmente. 
	TensorImageType::Pointer tensorImage=TensorImageType::New();
	tensorImage=tensorImageReader->GetOutput();
	
	if(m_VectorTensorData.size()>0){
		tensorImage->SetOrigin(m_VectorTensorData[0].image->GetOrigin());
	}
	this->DistribuirTensor( tensorImage, m_activeinput, m_filename );
	//this->DistribuirTensor( tensorImageReader->GetOutput(), m_activeinput, m_filename );
	////////////////////////////////////////////////////////////////////////////////////////////
	
	tensorImageReader = NULL;

	MyTensorConsole->SetHasBeenRegistered(false);
	MyTensorConsole->SetHasBeenFiberTracked(false);
}


void UsimagToolConsole::LoadStrainTensor(  )
{

	typedef itk::StrainTensorToSymTensorImageFilter<VectorImage4DType, STImageType>		CastFilterType;
	typedef CastFilterType::Pointer								CastFilterPointerType;
	
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}

	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	STReaderType::Pointer tensorImageReader = STReaderType::New();
	tensorImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	
	tensorImageReader->AddObserver( itk::ProgressEvent(), callback );

	// Se generan los nombres sustituyendo el "00" por un número de la serie
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	std::string format (m_filename);
	format.replace(format.find("00"), 2, "%02d");
	nameGenerator->SetSeriesFormat(format);

	nameGenerator->SetStartIndex(0);
	nameGenerator->SetEndIndex(100);
	nameGenerator->SetIncrementIndex(1);

	// Se obtiene el número de ficheros de la serie
	int numFiles = 0;
	while ( itksys::SystemTools::FileExists( nameGenerator->GetFileNames()[numFiles].c_str() ) ) {
		numFiles++;
	}

	// Se abren los ficheros de tensor de esfuerzo en primer lugar
	nameGenerator->SetStartIndex(0);
	nameGenerator->SetEndIndex((numFiles-1)/2);
	tensorImageReader->SetFileNames(nameGenerator->GetFileNames());

	try{
		tensorImageReader->Update();
	}
	catch( itk::ExceptionObject & e ) {
		fl_alert( e.GetDescription() );
		return;
	}

	STImageType *tensorImage = STImageType::New();

	// Se convierten los píxeles de tipo Vector a StrainTensor
	CastFilterPointerType caster = CastFilterType::New();
	caster->SetInput(tensorImageReader->GetOutput());
	tensorImage = caster->GetOutput();


	// En segundo lugar se leen los datos de deformación
	DeformReaderType::Pointer deformImageReader = DeformReaderType::New();
	deformImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	
	deformImageReader->AddObserver( itk::ProgressEvent(), callback );

	nameGenerator->SetStartIndex((numFiles-1)/2 + 1);
	nameGenerator->SetEndIndex(numFiles-1);
	nameGenerator->SetIncrementIndex(1);
	
	deformImageReader->SetFileNames(nameGenerator->GetFileNames());

	
	try{
		deformImageReader->Update();
	}
	catch( itk::ExceptionObject & e ) {
		fl_alert( e.GetDescription() );
		return;
	}

	DeformImageType *deformImage = DeformImageType::New();
	deformImage = deformImageReader->GetOutput();


	// Si cargamos dos volumenes se ajusta el origen para que tenga el mismo valor en las dos imágenes
	// y así poder compararlos mejor visualmente. 
	if(m_VectorSTData.size()>0){
		tensorImage->SetOrigin(m_VectorSTData[0].image->GetOrigin());
		tensorImage->SetSpacing(m_VectorSTData[0].image->GetSpacing());
	 	deformImage->SetOrigin(m_VectorSTData[0].image->GetOrigin());
		deformImage->SetSpacing(m_VectorSTData[0].image->GetSpacing());
	}

	tensorImage->Update();

	// El nombre se obtiene del nombre del fichero
	std::string nombre(m_filename);
	int pos = nombre.rfind('/');
	nombre = nombre.substr(nombre.rfind('/') + 1);

	// Se añade la nueva imagen al vector
	DataSTElementType tensor_data;
	int dataId = m_VectorSTData.size();
	char bname[200];
	tensor_data.Id     = dataId;
	tensor_data.nombre = nombre;
	tensor_data.image  = tensorImage;
	tensor_data.deform_image = deformImage;

	m_VectorSTData.push_back(tensor_data);

	tensorImageReader = NULL;

	// Se cambia la vista para mostrar la interfaz de tensor de esfuerzo
	this->ViewModeStrain3D();

}


/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadTIFF(  )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	FileReaderType::Pointer TIFFImageReader =  FileReaderType::New();
	TIFFImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	TIFFImageIOType::Pointer TIFFReader = TIFFImageIOType::New(); 
	TIFFImageReader->SetImageIO( TIFFReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	TIFFImageReader->SetFileName( m_filename );
	TIFFImageReader->AddObserver( itk::ProgressEvent(), callback );
	
	// Attempt to read
	try {
		TIFFImageReader->Update();
	}
	catch( itk::ExceptionObject & e ) {
		fl_alert( e.GetDescription() );
		return;
	}
	
	this->Distribuir( TIFFImageReader->GetOutput(), m_activeinput, m_filename );
	TIFFImageReader = NULL;
}


/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadPNG(  ) {
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	FileReaderType::Pointer PNGImageReader =  FileReaderType::New();
	PNGImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	PNGImageIOType::Pointer PNGReader = PNGImageIOType::New(); 
	PNGImageReader->SetImageIO( PNGReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	PNGImageReader->SetFileName( m_filename );
	PNGImageReader->AddObserver( itk::ProgressEvent(), callback );
	
	// Attempt to read
	try{
		PNGImageReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	this->Distribuir( PNGImageReader->GetOutput(), m_activeinput, m_filename );
	PNGImageReader = NULL;
}


/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadDicom(  )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	FileFloatReaderType::Pointer DicomImageReader = FileFloatReaderType::New();
	DicomImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	DicomImageIOType::Pointer DicomReader = DicomImageIOType::New(); 
	DicomImageReader->SetImageIO( DicomReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	DicomImageReader->SetFileName( m_filename );
	DicomImageReader->AddObserver( itk::ProgressEvent(), callback );
	
	try {
		DicomImageReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	this->Distribuir( DicomImageReader->GetOutput(), m_activeinput, m_filename );
	DicomImageReader = NULL;
}

void UsimagToolConsole::LoadGenericDWI(  )
{
	if( !m_filename ){return;}
	int l;
	for( l=strlen(m_filename)-1; l>0; --l ){
		if( m_filename[l]=='.' ){break;}
	}
	char tail1[6] = ".nhdr";
	char tail2[6] = ".nrrd";
	bool isNrrd = true;
	unsigned int j = l;
	for( int k=0; k<5; ++k, ++j ){
		if( j==strlen(m_filename) ){
			isNrrd = false;
			break;
		}
		if( tail1[k]!=m_filename[j] ){
			isNrrd = false;
			break;
		}
	}
	if( !isNrrd ){
		isNrrd = true;
		j = l;
		for( int k=0; k<5; ++k, ++j ){
			if( j==strlen(m_filename) ){
				isNrrd = false;
				break;
			}
			if( tail2[k]!=m_filename[j] ){
				isNrrd = false;
				break;
			}
		}
	}
	if( isNrrd )
		this->LoadNrrdDWI();
	else
		this->LoadDicomDWI();
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadGenericDWIandProcess(  )
{
	if( !m_filename ){return;}
	int l;
	for( l=strlen(m_filename)-1; l>0; --l ){
		if( m_filename[l]=='.' ){break;}
	}
	char tail1[6] = ".nhdr";
	char tail2[6] = ".nrrd";
	bool isNrrd = true;
	unsigned int j = l;
	for( int k=0; k<5; ++k, ++j ){
		if( j==strlen(m_filename) ){
			isNrrd = false;
			break;
		}
		if( tail1[k]!=m_filename[j] ){
			isNrrd = false;
			break;
		}
	}
	if( !isNrrd ){
		isNrrd = true;
		j = l;
		for( int k=0; k<5; ++k, ++j ){
			if( j==strlen(m_filename) ){
				isNrrd = false;
				break;
			}
			if( tail2[k]!=m_filename[j] ){
				isNrrd = false;
				break;
			}
		}
	}
	if( isNrrd )
		this->LoadNrrdDWIandProcess();
	else
		this->LoadDicomDWIandProcess();
}

void UsimagToolConsole::LoadNrrdDWI(  ) {
  std::cout << "Loading Nrrd DWI" << std::endl;	
  if( !m_filename ){return;}
  CommandPointer callback = CommandType::New();
  callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );

  NrrdToDWImagesReaderType::Pointer NrrdToDWImagesReader = NrrdToDWImagesReaderType::New();
  NrrdToDWImagesReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
  NrrdToDWImagesReader->SetFileName( m_filename );
  NrrdToDWImagesReader->AddObserver( itk::ProgressEvent(), callback );

  try{
	NrrdToDWImagesReader->Update();
  }
  catch( itk::ExceptionObject & e ){
	fl_alert( e.GetDescription() );
	return;
  }

  // Leemos los DWIs
  DWImagesType::Pointer DWImageData;
  DWImageData = NrrdToDWImagesReader->GetOutput();
  NrrdToDWImagesReader = NULL;

	std::string           nombre = "DWI";
	DataDWIElementType dwi_data;
	dwi_data.Id     = m_VectorDWIData.size();
	dwi_data.nombre = nombre.c_str();
	dwi_data.image  = DWImageData;
	m_VectorDWIData.push_back(dwi_data);
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadNrrdDWIandProcess(  )
{
	std::cout << "Leyendo y procesando DWI Nerrd" << std::cout;
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	NrrdToDWImagesReaderType::Pointer NrrdToDWImagesReader = NrrdToDWImagesReaderType::New();
	NrrdToDWImagesReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	NrrdToDWImagesReader->SetFileName( m_filename );
	NrrdToDWImagesReader->AddObserver( itk::ProgressEvent(), callback );
	
	try{
		NrrdToDWImagesReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	// Leemos los DWIs
	DWImagesType::Pointer DWImageData;
	DWImageData = NrrdToDWImagesReader->GetOutput();
	NrrdToDWImagesReader = NULL;
	
	// Filtrar los datos obtenidos
	try{
		MyTensorConsole->FilterDWI(DWImageData);
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	// Creamos espacio para el tensor
	TensorImageType::Pointer tensor;
	// Creamos espacio para la imagen t2
	InputImageType::Pointer  t2image;
	// Estimamos el tensor:
	try{
		MyTensorConsole->EstimateTensor( DWImageData, tensor, t2image, true );
		DWImageData = NULL;
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	// Metemos la imagen de tensores calculada en el tensordataBrowser:
	DataTensorElementType tensor_data;
	char bname[200];
	str_nombre( m_filename, bname );
	std::string nombre( bname );
	tensor_data.Id     = m_VectorTensorData.size();
	tensor_data.nombre = nombre.c_str();
	tensor_data.image  = tensor;
	
	m_VectorTensorData.push_back(tensor_data);
	
	// Metemos la imagen T2 calculada en el dataBrowser:
	DataElementType data_t2;
	std::string     nombre_t2 = "T2 " + nombre;
	data_t2.Id     = m_VectorData.size();
	data_t2.nombre = nombre_t2.c_str();
	data_t2.image  = t2image;
	
	m_VectorData.push_back(data_t2);
	// mostramos la T2 en el segundo viewer:
	ImageViewer[1]->SetImage(data_t2.image);
	ShowImageSrc(1);  
	
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	// Move this definition to UsimagToolBase.h!
	typedef itk::ComputeFAFilter<TensorImageType,InputImageType> FAfilterType;
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	FAfilterType::Pointer FAfilter = FAfilterType::New();
	FAfilter->SetInput(tensor);
	try{
		FAfilter->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	// Declaramos una imagen escalar donde poner la FA y la ponemos en el dataBrowser:
	DataElementType data;
	data.Id     = m_VectorData.size();
	data.nombre = "FA " + nombre;
	data.image  = FAfilter->GetOutput();
	
	m_VectorData.push_back(data);
	// Mostramos la imagen de FA en el segundo 
	ImageViewer[0]->SetImage(data.image);
	ShowImageSrc(0);  	
}

void UsimagToolConsole::LoadDicomDWI(  )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	char m_directory[300];
	directorio( m_filename, m_directory );
	
	DICOMtoDWImagesReaderType::Pointer DICOMtoDWImagesReader = DICOMtoDWImagesReaderType::New();
	DICOMtoDWImagesReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	DICOMtoDWImagesReader->SetDirectoryName( m_directory );
	DICOMtoDWImagesReader->AddObserver( itk::ProgressEvent(), callback );
	
	try{
		DICOMtoDWImagesReader->Update();
	} 
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}

	DWImagesType::Pointer DWImageData = DICOMtoDWImagesReader->GetOutput();
	DICOMtoDWImagesReader = NULL;
	std::string           nombre = "DWI";
	DataDWIElementType dwi_data;
	dwi_data.Id     = m_VectorDWIData.size();
	dwi_data.nombre = nombre.c_str();
	dwi_data.image  = DWImageData;
	m_VectorDWIData.push_back(dwi_data);
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadDicomDWIandProcess(  )
{
	if( !m_filename ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	char m_directory[300];
	directorio( m_filename, m_directory );
	
	DICOMtoDWImagesReaderType::Pointer DICOMtoDWImagesReader = DICOMtoDWImagesReaderType::New();
	DICOMtoDWImagesReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	DICOMtoDWImagesReader->SetDirectoryName( m_directory );
	DICOMtoDWImagesReader->AddObserver( itk::ProgressEvent(), callback );
	
	try{
		DICOMtoDWImagesReader->Update();
	} 
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** Is it necessary to update the filter twice? */
	/*DICOMtoDWImagesReader->Modified(); // This is only necessary if we need to update the filter twice
	// Leemos la serie Dicom de DWIs
	try{
		DICOMtoDWImagesReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}*/
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	//std::string           nombre      = DICOMtoDWImagesReader->GetPatientName();
	std::string           nombre = "Tensor";
	DWImagesType::Pointer DWImageData = DICOMtoDWImagesReader->GetOutput();
	DICOMtoDWImagesReader = NULL;
	
	// Filtrar los datos obtenidos
	try{
		MyTensorConsole->FilterDWI( DWImageData );
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	// Creamos espacio para el tensor
	TensorImageType::Pointer tensor;
	// Creamos espacio para la imagen t2
	InputImageType::Pointer  t2image;
	// Estimamos el tensor:
	try{
		MyTensorConsole->EstimateTensor( DWImageData, tensor, t2image, true );
		DWImageData = NULL;
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	// Metemos la imagen de tensores calculada en el tensordataBrowser:
	DataTensorElementType tensor_data;
	tensor_data.Id     = m_VectorTensorData.size();
	tensor_data.nombre = nombre.c_str();
	tensor_data.image  = tensor;
	
	m_VectorTensorData.push_back(tensor_data);
	/** This code is no longer necessary:
	m_tensordataBrowser->add(nombre.c_str());
	m_tensordataBrowser->select(m_tensordataBrowser->size());
	m_tensordataBrowser->redraw();
	*/
	
	// Metemos la imagen T2 calculada en el dataBrowser:
	DataElementType data_t2;
	std::string     nombre_t2;
	data_t2.Id     = m_VectorData.size();
	nombre_t2      = "T2 " + nombre;
	data_t2.nombre = nombre_t2.c_str();
	data_t2.image  = t2image;
	
	m_VectorData.push_back(data_t2);
	/** This code is no longer necessary:
	m_dataBrowser->add(nombre_t2.c_str());
	m_dataBrowser->select(m_dataBrowser->size());
	m_dataBrowser->redraw();
	*/
	
	/** This code is no longer necessary:
	// Actualizamos los Box Lists:
	m_destino->add(nombre_t2.c_str(),0,NULL,NULL,0);
	m_op1->add(nombre_t2.c_str(),0,NULL,NULL,0);
	m_op2->add(nombre_t2.c_str(),0,NULL,NULL,0);
	 */
	
	// mostramos la T2 en el segundo viewer:
	ImageViewer[1]->SetImage(data_t2.image);
	ShowImageSrc(1);
	
	
	
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	// Move this definition to UsimagToolBase.h!
	typedef itk::ComputeFAFilter<TensorImageType,InputImageType> FAfilterType;
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	FAfilterType::Pointer FAfilter = FAfilterType::New();
	FAfilter->SetInput(tensor);
	try{
		FAfilter->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	// Declaramos una imagen escalar donde poner la FA y la ponemos en el dataBrowser:
	DataElementType data;
	data.Id = m_VectorData.size();
	data.nombre = "FA " + nombre;
	data.image = FAfilter->GetOutput();
	
	m_VectorData.push_back(data);
	
	/** This code is no longer necessary:
	m_dataBrowser->add(nombre.c_str());
	m_dataBrowser->select(m_dataBrowser->size());
	m_dataBrowser->redraw();
	*/
	
	/** This code is no longer necessary:
	// Actualizamos los Box Lists:
	m_destino->add(nombre.c_str(),0,NULL,NULL,0);
	m_op1->add(nombre.c_str(),0,NULL,NULL,0);
	m_op2->add(nombre.c_str(),0,NULL,NULL,0);
	 */
	
	// Mostramos la imagen de FA en el segundo
	ImageViewer[0]->SetImage(data.image);
	ShowImageSrc(0);
}



/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadDicomSeries( )
{
	if( !m_filename ){return;}
    CommandPointer callback = CommandType::New();  
    callback->SetCallbackFunction(this,&UsimagToolConsole::OnProgress);

    // Falta poner: m_directory, y m_series_name (opcional)
    char m_directory[300];
    directorio( m_filename, m_directory );
    //char *m_series_name;
    int series_name = 0;  
    
    std::cout << std::endl << "The directory: " << std::endl;
    std::cout << std::endl << m_directory << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;

	itk::DICOMSeriesFileNames::Pointer nameGenerator = itk::DICOMSeriesFileNames::New();
    nameGenerator->SetDirectory( m_directory );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	// Move this definition to UsimagToolBase.h!
    typedef std::vector<std::string> seriesIdContainer;
    /** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	const seriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

    seriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    seriesIdContainer::const_iterator seriesEnd = seriesUID.end();
  
    while( seriesItr != seriesEnd ) {
		std::cout << seriesItr->c_str() << std::endl;
		seriesItr++;
    }
  
    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;

	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	// Move this definition to UsimagToolBase.h!
    typedef std::vector<std::string> fileNamesContainer;
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
    fileNamesContainer fileNames;

    if( series_name == 0 ) { // If no optional third argument 
		std::cout << seriesUID.begin()->c_str() << std::endl;
		fileNames = nameGenerator->GetFileNames();
    }
	else{
		//fileNames = m_nameGenerator->GetFileNames( m_series_name );
		//std::cout << m_series_name << std::endl;
	}
    std::cout << std::endl << std::endl;

    SeriesReaderType::Pointer DicomSeriesImageReader = SeriesReaderType::New();
	DicomSeriesImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	DicomSeriesImageIOType::Pointer DicomSeriesReader = DicomSeriesImageIOType::New();
	DicomSeriesImageReader->SetImageIO( DicomSeriesReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	DicomSeriesImageReader->SetFileNames( fileNames );    
    DicomSeriesImageReader->AddObserver( itk::ProgressEvent(), callback );
    
	try{
		DicomSeriesImageReader->Update();
    }
	catch (itk::ExceptionObject &ex) {
		fl_alert( ex.GetDescription() );
		return;
    }

    this->Distribuir( DicomSeriesImageReader->GetOutput(), m_activeinput, m_filename );
	DicomSeriesImageReader = NULL;
}


/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::LoadNrrd(  ) {
	if( !m_filename ){return;}
    CommandPointer callback = CommandType::New();  
    callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	
	FileReaderType::Pointer NrrdImageReader = FileReaderType::New();
	NrrdImageReader->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** TO DO: Check if this piece of code is really necessary. */
	NrrdImageIOType::Pointer NrrdReader = NrrdImageIOType::New();
	NrrdImageReader->SetImageIO( NrrdReader );
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	NrrdImageReader->SetFileName( m_filename );
	NrrdImageReader->AddObserver( itk::ProgressEvent(), callback );
	
	// Attempt to read
	try{
		NrrdImageReader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	this->Distribuir( NrrdImageReader->GetOutput(), m_activeinput, m_filename );
	NrrdImageReader = NULL;
}



void UsimagToolConsole::LoadModel() {
  m_filename = fl_file_chooser("Image Filename","*.{vtk}","");
  if( !m_filename ){return;}
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
  reader->SetFileName(m_filename);
  try{
	  reader->Update();
  }
  catch( itk::ExceptionObject & e ){
	  fl_alert( e.GetDescription() );
	  return;
  }
  vtkPolyData* data = vtkPolyData::New();
  data->ShallowCopy( reader->GetOutput() );
  vtkActor *Actor   = vtkActor::New();
  DataModelElementType model_data;
  model_data.Id     = m_VectorModelData.size();
  char nombre[200];
  str_nombre( m_filename, nombre );
  model_data.nombre = nombre;
  model_data.data   = data;
  model_data.actor  = Actor;
  m_VectorModelData.push_back( model_data );
  ImageViewer3D->ConnectMapperFiber( data, Actor );
}

void UsimagToolConsole::LoadFibers(void) {
	m_filename = fl_file_chooser("Image Filename","*.{vtk}","");
	if( !m_filename ){return;}
	vtkPolyDataReader* reader=vtkPolyDataReader::New();
	reader->SetFileName(m_filename);
	try{
		reader->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	vtkPolyData* fibers = vtkPolyData::New();
	
	vtkActor *Actor = vtkActor::New();
	DataModelElementType model_data;
    model_data.Id = m_VectorModelData.size();
    model_data.nombre = "Streamline";
    model_data.data = fibers;
    model_data.actor = Actor;
    m_VectorModelData.push_back( model_data );
	double color[3] = {1,1,1};
	ImageViewer3D->renderTracts( fibers, Actor, color, 0.3 );
    std::cout << "rendering fibers" << std::endl;
	//reader->Delete();
} 

void UsimagToolConsole::WriteFibers(vtkPolyData* streamlines) {
    m_filename = fl_file_chooser("Image Filename","*.{vtk}","");
	if( !m_filename ){return;}
	vtkPolyDataWriter* writer=vtkPolyDataWriter::New();
	writer->SetFileTypeToASCII();
	writer->SetInput(streamlines);	
	writer->SetFileName(m_filename);
	try{
		writer->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
}

void UsimagToolConsole::DistribuirTensor( TensorImageType::Pointer aux, int window, const char* filename)
{
	char nombre[200];
	str_nombre( filename, nombre );
	char nombre_tensor[200];
	sprintf( nombre_tensor, "TENSOR_%s", nombre );
	
	DataTensorElementType tensor_data;
	tensor_data.Id     = m_VectorTensorData.size();
	tensor_data.nombre = nombre_tensor;
	tensor_data.image  = aux;
	
	m_VectorTensorData.push_back(tensor_data);
	/** This code is no longer necessary:
	m_tensordataBrowser->add(nombre_tensor);
	m_tensordataBrowser->select(m_tensordataBrowser->size());
	m_tensordataBrowser->redraw();
	*/
	
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	// TODO: Move this typedef!!!
	typedef itk::ComputeFAFilter<TensorImageType,InputImageType> FAfilterType;
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	FAfilterType::Pointer FAfilter = FAfilterType::New();
	FAfilter->SetInput(aux);
	try{
		FAfilter->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}

	sprintf(nombre,"%s_FA",nombre);
	DataElementType data;
	data.Id      = m_VectorData.size();
	data.nombre  = nombre;
	data.image   = FAfilter->GetOutput();
	m_VectorData.push_back(data);
	
	/** This code is no longer necessary:
	m_dataBrowser->add(nombre);
	m_dataBrowser->select(m_dataBrowser->size());
	m_dataBrowser->redraw();
	*/
	/** This code is no longer necessary:
	m_destino->add(nombre,0,NULL,NULL,0);
	m_op1->add(nombre,0,NULL,NULL,0);
	m_op2->add(nombre,0,NULL,NULL,0);
	*/
	
	std::cout << aux->GetDirection() << std::endl;
	
	ImageViewer[window]->SetImage( data.image );
	float v[3][3];
	for (unsigned int i=0;i<3;i++) {
	  for (unsigned int j=0;j<3;j++) {
	    v[i][j] = (float)aux->GetDirection()[i][j];
	  }
	}
	ImageViewer[window]->SetOrientationMatrix( &v[0][0] );
	ShowImageSrc(m_activeinput);
}


/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::Distribuir( InputImageType::Pointer aux, int window, const char* filename)
{
	char nombre[200];
	str_nombre( filename, nombre );
	std::cout << aux->GetDirection() << std::endl;
	m_orientation = check_orientation(aux->GetDirection() );
	//OrientarImagen(aux, m_orientation);
	DataElementType data;
	data.Id     = m_VectorData.size();
	data.nombre = nombre;
	data.image  = aux;
	m_VectorData.push_back( data );
	
	/** This code is no longer necessary:
	m_dataBrowser->add(nombre);  
	m_dataBrowser->select(m_dataBrowser->size());
	m_dataBrowser->redraw();
	*/
	
	/* This code is no longer necessary:
	m_destino->add( nombre, 0, NULL, NULL, 0 );
	m_op1->add(     nombre, 0, NULL, NULL, 0 );
	m_op2->add(     nombre, 0, NULL, NULL, 0 );
	*/
	float v[3][3];
	for (unsigned int i=0;i<3;i++) {
	  for (unsigned int j=0;j<3;j++) {
	    v[i][j] = (float)aux->GetDirection()[i][j];
	  }
	}
	
	ImageViewer[window]->SetImage( data.image );
	ImageViewer[window]->SetOrientationMatrix( &v[0][0] );

	ShowImageSrc( m_activeinput );  
}

int UsimagToolConsole::check_orientation(InputImageType::DirectionType Dir)
{

    if( (abs(Dir[0][2]) > abs(Dir[1][2])) && (abs(Dir[0][2]) > abs(Dir[2][2]))) {
		return 0;
	}
	if( (abs(Dir[1][2]) > abs(Dir[0][2])) && (abs(Dir[1][2]) > abs(Dir[2][2]))) {
		return 1;
 	}
    if( (abs(Dir[2][2]) > abs(Dir[0][2])) && (abs(Dir[2][2]) > abs(Dir[1][2]))) {
		return 2;
	}
	return -1;
}

void UsimagToolConsole::OrientarImagen( InputImageType::Pointer &image, int orientation )
{
  typedef itk::OrientImageFilter<InputImageType,InputImageType> OrientImageFilterType; 
  OrientImageFilterType::Pointer orienter = OrientImageFilterType::New(); 

  orienter->SetInput( image );
  orienter->SetUseImageDirection( true );
  switch (orientation){
    case 0:
	  orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL); 
	  break;
    case 1:
	  orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP); 
	  break;
	case 2:
	  orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI); 
      break;
	default:
	  break;
  }
  try {
	  orienter->Update();
	  orienter->Modified();
  } 
  catch(itk::ExceptionObject &excp){
	  fl_alert( excp.GetDescription() );
	  return;
  } 
	
  image = orienter->GetOutput();
 
  return; 
} 


void UsimagToolConsole::DistribuirColor( ImageRGBType::Pointer aux, int window, const char* filename )
{
	char nombre[200];
	str_nombre( filename, nombre );

	
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	// It makes no sense to add an entry to m_VectorData, since it is not the same kind as the data in m_VectorData!!!
	m_NImagesCargadas++;
	/*
	DataElementType data;
	data.Id     = m_VectorData.size();
	data.nombre = nombre;
	data.image  = aux;
	m_VectorData.push_back(data);
	*/
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	
	/** This code is no longer necessary:
	m_dataBrowser->add(nombre);  
	m_dataBrowser->select( m_dataBrowser->size() );
	m_dataBrowser->redraw();
	*/
	
	/** This code is no longer necessary:
	m_destino->add( nombre, 0, NULL, NULL, 0 );
	m_op1->add(     nombre, 0, NULL, NULL, 0 );
	m_op2->add(     nombre, 0, NULL, NULL, 0 );
	*/
	
	ImageColorViewer->SetInputImage( aux );
	ShowImageColor(3);  
}

void UsimagToolConsole::ShowImageColor(int window)
{
  ImageColorViewer->Show();
  ImageColorViewer->ViewDetails(false);
  ImageColorViewer->ViewValue(true);

  /////////////////
  
  sliceNumberSlider[window]->range( 0.0f, ImageColorViewer->numSlices()-1 );
  sliceNumberSlider[window]->value((float)ImageColorViewer->sliceNum());
  
  ImageColorViewer->update();
  glWindowGroup[window]->redraw();
  Botones[window]->redraw();
  Fl::check();
  ImageColorViewer->redraw();
  Fl::check();
}

void UsimagToolConsole::ShowImageSrc(int window)
{
  ImageViewer[window]->Show();
  ImageViewer[window]->ViewDetails(false);
  ImageViewer[window]->ViewValue(false); 
  unsigned int midvalue = (unsigned int)(ImageViewer[window]->numSlices()-1)/2;
  ImageViewer[window]->SelectSlice(midvalue);

  switch (m_orientation) {
    case 0:
         ImageViewer[window]->SetSagital();
		 break;
	case 1:
	     ImageViewer[window]->SetCoronal();
		 break;
	case 2:
	     ImageViewer[window]->SetAxial();
		 break;
	default:
		 break;
  }
  
  ///// Para los controles de intensidad de imagen: 
  // Mi metodo:
  float dataMax = ImageViewer[window]->dataMax();
  float dataMin = ImageViewer[window]->dataMin();
  float iwDiff  = dataMax - dataMin;
  double iwStep = (dataMax-dataMin)/50.0;
  sliderIW[window]->range(dataMin-iwDiff*2,dataMin+iwDiff*2);
  sliderIL[window]->range(dataMax-iwDiff*2,dataMax+iwDiff*2);
	
  sliderIW[window]->step(iwStep);
  sliderIL[window]->step(iwStep);
  sliderIW[window]->value(dataMin);
  sliderIL[window]->value(dataMax); 
  
  m_minintensity_output[window]->value(dataMin);
  m_maxintensity_output[window]->value(dataMax);
  /////////////////
	sliderCMin[window]->range(dataMin,dataMax);
	sliderCMin[window]->value(dataMin);
	sliderCMin[window]->step((dataMax-dataMin)/100);
	
	sliderCMax[window]->range(dataMin,dataMax);
	sliderCMax[window]->value(dataMax);
	sliderCMax[window]->step((dataMax-dataMin)/100);
	
	
  sliceNumberSlider[window]->range( 0.0f, ImageViewer[window]->numSlices()-1 );
  sliceNumberSlider[window]->value((float)ImageViewer[window]->sliceNum());
   
  ImageViewer[window]->update();
  glWindowGroup[window]->redraw();
  Botones[window]->redraw();
  sliderIW[window]->redraw();
  sliderIL[window]->redraw();
	//////////////////////////
	//////////////////////////
	//////////////////////////
    // Revisar esto!
	InitializeOverlay(window);
	//////////////////////////
	//////////////////////////
	//////////////////////////
  ///////////////////////////////
  Fl::check();
  ImageViewer[window]->redraw();
  Fl::check();
}

void UsimagToolConsole::ViewData(int window, int dataId) {
  if (window == 3) return;
  if (m_VectorData[dataId].image != static_cast<InputImageType::Pointer>(NULL) ) {
    m_orientation = check_orientation(m_VectorData[dataId].image->GetDirection() );
	  
    ImageViewer[window]->SetImage( m_VectorData[dataId].image);
	this->ShowImageSrc(window);
	float v[3][3];
	  for (unsigned int i=0;i<3;i++) {
		  for (unsigned int j=0;j<3;j++) {
	  	    v[i][j] = (float)m_VectorData[dataId].image->GetDirection()[i][j];
	  	  }
	  }
	ImageViewer[window]->SetOrientationMatrix( &v[0][0] );
   }
}

void UsimagToolConsole::LoadOverlay(int window) {
	m_filename = fl_file_chooser("Image Filename","*","");
	if( !m_filename ){return;}
	FileUCharReaderType::Pointer reader = FileUCharReaderType::New();
	reader->SetFileName(m_filename);
	reader->Update();
	typedef itk::ImageRegionConstIterator<UCharImageType> UCharConstIteratorType;
	typedef itk::ImageRegionIterator<UCharImageType> UCharIteratorType;
	
	if (ImageViewer[window]->GetInputImage() != static_cast<InputImageType::Pointer>(NULL) ) {
		if (ImageViewer[window]->GetInputImage()->GetLargestPossibleRegion().GetSize()  == reader->GetOutput()->GetLargestPossibleRegion().GetSize() ) {
		   m_ImageOverlay->SetRegions( ImageViewer[window]->GetInputImage()->GetLargestPossibleRegion() );
	       m_ImageOverlay->Allocate();
			UCharIteratorType it1(m_ImageOverlay, m_ImageOverlay->GetRequestedRegion());
			UCharConstIteratorType it2(reader->GetOutput(), reader->GetOutput()->GetRequestedRegion());
			//m_ImageOverlay = reader->GetOutput();
			for (it1.GoToBegin(), it2.GoToBegin();!it1.IsAtEnd();++it1,++it2) {
				it1.Set(it2.Get());
			}
		}
	}
	SetOverlay(window);
}


void UsimagToolConsole::InitializeOverlay(int window) {
	if (ImageViewer[window]->GetInputImage() != static_cast<InputImageType::Pointer>(NULL) ) {
	    m_ImageOverlay->SetRegions( ImageViewer[window]->GetInputImage()->GetLargestPossibleRegion() );
	    m_ImageOverlay->Allocate();
	    m_ImageOverlay->FillBuffer(0);
	}
	SetOverlay(window);
}

void UsimagToolConsole::SetOverlay(int window) {
   if (ImageViewer[window]->GetInputImage()->GetLargestPossibleRegion().GetSize()  == m_ImageOverlay->GetLargestPossibleRegion().GetSize() ) {
     ImageViewer[window]->SetOverlay( m_ImageOverlay );
     ImageViewer[window]->SetOverlayOpacity(0.8);
     ImageViewer[window]->ViewOverlayData(true);
     ImageViewer[window]->Update();
   }
}
	

void UsimagToolConsole::ClearOverlay( ) {
	if (m_ImageOverlay != static_cast<UCharImageType::Pointer>(NULL) ) {
	  m_ImageOverlay->FillBuffer(0);		
	  for (unsigned int i=0;i<3;i++) 
	    ImageViewer[i]->Update();
	}
}


	
void UsimagToolConsole::ViewSliceIn3DWindow(InputImageType::Pointer image, int orientation) {
  vtkImageData* data;  
  vtkImageImport* vtkImporter = vtkImageImport::New();
  m_VTKexporter->SetInput(image);
  ConnectPipelines(m_VTKexporter, vtkImporter);

  try{
	  vtkImporter->GetOutput()->Update();
  }
  catch(itk::ExceptionObject &excp){
	  fl_alert( excp.GetDescription() );
	  return;
  }
  
  data = vtkImporter->GetOutput();
  
  ImageViewer3D->renderSlice( data, orientation );
  
  Fl::check();
  ImageViewer3D->redraw();
  Fl::check();
}


void UsimagToolConsole::ViewStrainSlice3D() {

/***************************************************************************************/
	
	int dataId = m_strainTensorDataBrowser->value()-1;
	if( dataId<0){return;}

	STImageType::Pointer image = (m_VectorSTData)[dataId].image;
	ExtractFilterType::Pointer extract = ExtractFilterType::New();
	extract->SetInput( image );
	STImageType::RegionType inputRegion = image->GetLargestPossibleRegion();
	STImageType::SizeType inputSize = inputRegion.GetSize();
	inputSize[3] = 0;
	STImageType::IndexType inputIndex = inputRegion.GetIndex();
	inputIndex[3] = strainTime->value(); // Instante temporal
	STImageType::RegionType desiredRegion;
	desiredRegion.SetSize( inputSize );
	desiredRegion.SetIndex( inputIndex );
	extract->SetExtractionRegion( desiredRegion );

	ComputeStrainScalarsType::Pointer filter = ComputeStrainScalarsType::New();
	filter->SetInput( extract->GetOutput() );

	filter->SetComputeInvariant();

	try{
		filter->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}

	char nombre[200];
	sprintf( nombre, "INV" );
	
  	vtkImageImport* vtkImporter = vtkImageImport::New();

  	m_VTKexporter->SetInput(filter->GetOutput());
  	ConnectPipelines(m_VTKexporter, vtkImporter);

  	try{
		  vtkImporter->GetOutput()->Update();
  	}
  	catch(itk::ExceptionObject &excp){
	  	fl_alert( excp.GetDescription() );
	  	return;
  	}

  	this->m_activeStrainImage = vtkImporter->GetOutput();

  	this->m_activeStrainImage->Update();

  	ImageViewerStrain3D->renderSlice( this->m_activeStrainImage, 2, strainSlice->value() );
  
  	Fl::check();
  	ImageViewerStrain3D->redraw();
  	Fl::check();
}
  
void UsimagToolConsole::DoModel(int value1, int value2, int smooth_iter, InputImageType::Pointer image) {
  vtkImageData* data;
  
  vtkImageImport* vtkImporter = vtkImageImport::New();
  m_VTKexporter->SetInput(image);
  ConnectPipelines(m_VTKexporter, vtkImporter);
  try{
	  vtkImporter->GetOutput()->Update();
  }
  catch(itk::ExceptionObject &excp){
	  fl_alert( excp.GetDescription() );
	  return;
  }
  data = vtkImporter->GetOutput();
  

  // Nuevo Pipeline:
  vtkContourFilter *modelExtractor = vtkContourFilter::New();
  modelExtractor->SetInput( data );
  modelExtractor->SetValue(value1, value2);
  vtkActor *Actor = vtkActor::New();
	
  vtkPolyData *newdata = vtkPolyData::New();

  vtkSmoothPolyDataFilter* smoother = vtkSmoothPolyDataFilter::New();
  smoother->SetInput(modelExtractor->GetOutput());
  smoother->BoundarySmoothingOn();
  smoother->SetRelaxationFactor((float)m_relaxfactor->value());
  smoother->SetNumberOfIterations((int)m_modeliter->value());
 
  vtkDecimatePro* decimate = vtkDecimatePro::New();
  decimate->SetInput(smoother->GetOutput());
  decimate->SetTargetReduction((float)m_targetreduction->value());
  vtkPolyDataNormals* normals = vtkPolyDataNormals::New();
  normals->SetInput(decimate->GetOutput());
  normals->FlipNormalsOn();

	try{
		normals->Update();
	}
	catch(itk::ExceptionObject &excp){
		fl_alert( excp.GetDescription() );
		return;
	}
		
  newdata->ShallowCopy(normals->GetOutput());
	
	
  DataModelElementType model_data;
  model_data.Id     = m_VectorModelData.size();
  model_data.nombre = "New";
  model_data.data   = newdata;
  model_data.actor  = Actor;
  m_VectorModelData.push_back( model_data );

  ImageViewer3D->ConnectMapper( newdata, Actor );
  std::cout << "Done" << std::endl;
	/** This code is no longer necessary:
  m_modeldataBrowser->add("Nuevo");
  m_modeldataBrowser->select(m_modeldataBrowser->size());
  m_modeldataBrowser->redraw();
	 */

}


void UsimagToolConsole::OnSliceChange(unsigned int value, unsigned int window) {
  // Anadido por la ventana a color
  if (window == 3) {
	ImageColorViewer->SelectSlice(value);
	sliceNumberSlider[3]->value(value);
    if (m_viewmode != 2) {
	  for (unsigned int i=0;i<3;i++) { 
        if ( ImageViewer[i]->GetInputImage()->GetLargestPossibleRegion().GetSize() == ImageColorViewer->GetInputImage()->GetLargestPossibleRegion().GetSize()) {
           ImageViewer[i]->SelectSlice(value);
           sliceNumberSlider[i]->value(value);
        }
	  }
	}
	return;
  }

  ImageViewer[window]->SelectSlice(value);
  if (m_viewmode != 2) {
    for (unsigned int i=0;i<3;i++) { 
	  if (i != window) {
        if ( ImageViewer[i]->GetInputImage()->GetLargestPossibleRegion().GetSize() == ImageViewer[window]->GetInputImage()->GetLargestPossibleRegion().GetSize() ){
          ImageViewer[i]->SelectSlice(value);
          sliceNumberSlider[i]->value(value);
        }
		if ( ImageColorViewer->GetInputImage()->GetLargestPossibleRegion().GetSize() == ImageViewer[window]->GetInputImage()->GetLargestPossibleRegion().GetSize() ){
			  ImageColorViewer->SelectSlice(value);
			  sliceNumberSlider[3]->value(value);
		}
      }
    }
  }
}

void UsimagToolConsole::OnOrientationChange(unsigned int window) {   
   if (window == 3) {
     ImageColorViewer->SetOrientation(m_Orientation[window]->value());
     sliceNumberSlider[window]->range( 0.0f, ImageColorViewer->numSlices()-1 );
	 return;
   }
   
   ImageViewer[window]->SetOrientation(m_Orientation[window]->value());
   sliceNumberSlider[window]->range( 0.0f, ImageViewer[window]->numSlices()-1 );
   unsigned int midvalue = (unsigned int)(ImageViewer[window]->numSlices()-1)/2;
   ImageViewer[window]->SelectSlice(midvalue);
   sliceNumberSlider[window]->value((float)ImageViewer[window]->sliceNum());
   if (m_viewmode != 2) {
   for (unsigned int i=0;i<3;i++) { 
	 if (i != window) {
       if ( ImageViewer[window]->GetInputImage()->GetLargestPossibleRegion().GetSize() == ImageViewer[i]->GetInputImage()->GetLargestPossibleRegion().GetSize() ){
	     ImageViewer[i]->SetOrientation(m_Orientation[window]->value());
	     sliceNumberSlider[i]->range( 0.0f, ImageViewer[window]->numSlices()-1 );
	     m_Orientation[i]->value(m_Orientation[window]->value());
		 ImageViewer[window]->SelectSlice(midvalue);
		 sliceNumberSlider[window]->value((float)ImageViewer[window]->sliceNum());
       }
	 }
   }
   }
}

void UsimagToolConsole::OnImageModeChange() {
  
  switch(m_ImageMode->value()){
		case 0: 
			ImageViewer[m_activeinput]->ImageMode(IMG_VAL);
			break;
		case 1: 
			ImageViewer[m_activeinput]->ImageMode(IMG_INV);
			break;
		case 2: 
			ImageViewer[m_activeinput]->ImageMode(IMG_LOG);
			break;
		case 3:	
			ImageViewer[m_activeinput]->ImageMode(IMG_DX);
			break;
		case 4: 
			ImageViewer[m_activeinput]->ImageMode(IMG_DY);
			break;
		case 5: 
			ImageViewer[m_activeinput]->ImageMode(IMG_DZ);
			break;
		case 6: 
			ImageViewer[m_activeinput]->ImageMode(IMG_BLEND);
			break;
		case 7:	
			ImageViewer[m_activeinput]->ImageMode(IMG_MIP);
			break;
		case 8: 
			ImageViewer[m_activeinput]->ImageMode(IMG_COLOR);
			break;
	}

}

void UsimagToolConsole::ChangeViewerColorMode() {
	ImageViewer[m_activeinput]->ChangeColorMode(m_colorMode->value());
	ImageViewer[m_activeinput]->update();
}


void UsimagToolConsole::OnCheckButtonChange(bool t, int window) {
  if (t) {
     m_activeinput = window;
	 for (int i=0;i<4;i++) {
	   if (i!=window) {
         m_checkbutton[i]->value(false);
       }
	 }
  }
}

void UsimagToolConsole::OnAmplify(Fl_Group* viewSplit, int window) {
  for (int i=0;i<3;i++) {
    if (i != window) {
      ImageViewer[i]->hide();
	  glWindowGroup[i]->hide();
      Botones[i]->hide();
	}
  }
  glWindowGroup[3]->hide();
  Botones[3]->hide();
  ImageColorViewer->hide();
  
  int x = viewSplit->x();
  int y = viewSplit->y();
  int w = viewSplit->w();
  int h = viewSplit->h();
  
  int wBotones = Botones[window]->w();
  int hBotones = Botones[window]->h();
  //Botones[window]->resize(426,26,400,20);
  Botones[window]->resize(x,y,wBotones,hBotones);
  m_Dism[window]->show();
  m_Ampl[window]->hide();
  //ImageViewer[0]->show();
  //glWindowGroup[window]->resize(410,446,800,650);
  glWindowGroup[window]->resize(x,y+hBotones,w,h-hBotones);
  ImageViewer[window]->resize(x,y+hBotones,w,h-hBotones-35);
  //sliceNumberSlider[window]->resize(x,y+h-hBotones-20,w,20);
  Fl::check();
  ImageViewer[0]->redraw();
  ImageViewer[1]->redraw();
  ImageViewer[2]->redraw();
  ImageColorViewer->redraw();
  //ImageViewer[3]->redraw();
  Fl::check();
}

void UsimagToolConsole::OnMinimize(int window,int x, int y) {
 for (int i=0;i<3;i++) {
    if (i != window) {
	  glWindowGroup[i]->show();
      ImageViewer[i]->show();
      Botones[i]->show();
    }
 }
 glWindowGroup[3]->show();
 ImageColorViewer->show();
 Botones[3]->show();
 Botones[window]->resize(x,y,400,20);
 glWindowGroup[window]->resize(x,y+1,400,340);
 m_Ampl[window]->show();
 m_Dism[window]->hide();
 Fl::check();
 ImageViewer[0]->redraw();
 ImageViewer[1]->redraw();
 ImageViewer[2]->redraw();
 Fl::check();
}

void UsimagToolConsole::UpdateIntensityW(int window, float val) {
  m_minintensity_output[window]->value(val);
  ImageViewer[window]->SetIntensityWindowingMin(val);
}

void UsimagToolConsole::UpdateIntensityL(int window, float val) {
  m_maxintensity_output[window]->value(val);
  ImageViewer[window]->SetIntensityWindowingMax(val);
}


void UsimagToolConsole::str_nombre(const char* file,char* nombre) {
  char aux[100];
  if (strrchr(file,'/') != NULL) {
    strcpy(aux,strrchr(file,'/'));
    strncpy(nombre,&aux[1],strlen(aux));
  }
  else{
    strcpy(nombre, file);
  }
}

void UsimagToolConsole::directorio(const char* file,char* dir) {
  int diferencia;
  diferencia = strlen(file);
  if (strrchr(file,'/') != NULL){
	  diferencia -= strlen(strrchr(file,'/'));
	  strncpy(dir,file,diferencia+1);
	  dir[diferencia] = '\0';
  }
  else{
    dir[0] = '.';
    dir[1] = '\0';
  }
}

void UsimagToolConsole::InfoImagen( ) {  
  int dataId = m_dataBrowser->value()-1;
  if (dataId < 0) return;
  char texto[100];
  float v[3][3];
  for (unsigned int i=0;i<3;i++) {
	  for (unsigned int j=0;j<3;j++) {
	    v[i][j] = (float)m_VectorData[dataId].image->GetDirection()[i][j];
	  }
  }
  
  Fl_Text_Buffer *buff = new Fl_Text_Buffer();
  sprintf(texto,"Orientation \n %f, %f, %f \n %f %f %f \n %f %f %f",v[0][0],v[0][1],v[0][2],v[1][0],v[1][1],v[1][2],v[2][0],v[2][1],v[2][2]);
  buff->text(texto);
  m_OrientationText->buffer(buff);
  
  char texto2[100];
  float origen[3];
  origen[0] = (float)m_VectorData[dataId].image->GetOrigin()[0];
  origen[1] = (float)m_VectorData[dataId].image->GetOrigin()[1];
  origen[2] = (float)m_VectorData[dataId].image->GetOrigin()[2];
  Fl_Text_Buffer *buff2 = new Fl_Text_Buffer();
  sprintf(texto2,"Origin %f, %f, %f",origen[0],origen[1],origen[2]);
  buff2->text(texto2);
  m_OriginText->buffer(buff2);
  
  char texto3[100];
  typedef itk::MinimumMaximumImageCalculator< InputImageType > MinMaxCalculatorType;
  MinMaxCalculatorType::Pointer calculator = MinMaxCalculatorType::New( );
  calculator->SetImage( m_VectorData[dataId].image );
  //InputImageType::PixelType minIntensity = calculator->GetMinimum( );
  //InputImageType::PixelType maxIntensity = calculator->GetMaximum( );
  calculator->Compute( );
  float minIntensity = calculator->GetMinimum( );
  float maxIntensity = calculator->GetMaximum( );
  std::cout << "maxIntensity " << maxIntensity << " minIntensity " << minIntensity << std::endl;  
  Fl_Text_Buffer *buff3 = new Fl_Text_Buffer();
  sprintf(texto3,"Min %f, Max %f",minIntensity,maxIntensity);
  buff3->text(texto3);
  m_MaxMinText->buffer(buff3);
 
  panelMessage->show();
}

/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::DeleteScalar( int dataId )
{
    if( dataId<(int)m_VectorData.size() ){
		//m_VectorData[dataId].image = NULL;
		m_VectorData.erase( dataId );
	}
}
void UsimagToolConsole::RenameScalar( int id, const char* name )
{
	if( id<0 || id>=(int)m_VectorData.size() )
		return;
	m_VectorData.rename( id, name );
}


/** Added: 11/07/2008 AN-TONIO */
void UsimagToolConsole::DeleteTensor( int dataId )
{
	if( dataId<(int)m_VectorTensorData.size() ){
		//m_VectorTensorData[dataId].image = NULL;
		m_VectorTensorData.erase( dataId );
	}
}

void UsimagToolConsole::RenameTensor( int id, const char* name )
{
	if( id<0 || id>=(int)m_VectorTensorData.size() )
		return;
	m_VectorTensorData.rename( id, name );
}



/** Last modified: 11/07/2008 AN-TONIO */
void UsimagToolConsole::DeleteModel( int dataId )
{
    std::cout << "borrando model " << dataId << std::endl;
	if( dataId<(int)m_VectorModelData.size() ){
		m_VectorModelData.erase( dataId );
	}
}

void UsimagToolConsole::RenameModel( int id, const char* name )
{
	if( id<0 || id>=(int)m_VectorModelData.size() )
		return;
	m_VectorModelData.rename( id, name );
}

void UsimagToolConsole::DeleteStrainTensor( int dataId )
{
	if( dataId<(int)m_VectorSTData.size() ){
		//m_VectorTensorData[dataId].image = NULL;
		m_VectorSTData.erase( dataId );
	}
}

void UsimagToolConsole::RenameStrainTensor( int id, const char* name )
{
	if( id<0 || id>=(int)m_VectorSTData.size() )
		return;
	m_VectorSTData.rename( id, name );
}

void UsimagToolConsole::ViewMode4x2D() {
  // variable m_viewmode para saber en que modo estamos
  m_viewmode = 0;
  m_view3D->hide();
  m_viewStrain3D->hide();
  m_viewSplit->show();
  ImageViewer[0]->Show();
  ImageViewer[1]->Show();
  ImageViewer[2]->Show();
  glWindowGroup[3]->show();
  //ImageViewer[3]->Show();
  ImageColorViewer->Show();
  Botones[3]->show();
  Fl::check();
  ImageViewer[0]->redraw();
  ImageViewer[1]->redraw();
  ImageViewer[2]->redraw();
  //ImageViewer[3]->redraw();
  ImageColorViewer->redraw();
  Fl::check();  
}

void UsimagToolConsole::ViewMode3D() {
  // variable m_viewmode para saber en que modo estamos
  m_viewmode = 1;
  m_viewSplit->hide();
  ImageViewer[0]->Hide();
  ImageViewer[1]->Hide();
  ImageViewer[2]->Hide();
  //ImageViewer[3]->Hide();
  ImageColorViewer->Hide();
  m_viewStrain3D->hide();
  ImageViewer3D->ResetView();
  m_view3D->resize(410,25,800,682);
  m_view3D->show();
  Fl::check();
  ImageViewer3D->redraw();
  Fl::check();
}

void UsimagToolConsole::ViewMode3_1() {
  // variable m_viewmode para saber en que modo estamos
  m_viewmode = 2;
  m_view3D->resize(810,366,400,320);
  m_viewSplit->show();
  ImageViewer[0]->Show();
  ImageViewer[1]->Show();
  ImageViewer[2]->Show();
  glWindowGroup[3]->hide();
  Botones[3]->hide();
  m_viewStrain3D->hide();
  m_view3D->show();
  /// Ponemos modo axial, sagital, coronal
  int Id = m_dataBrowser->value()-1;
  if (Id >= 0) {
    ImageViewer[0]->SetImage(m_VectorData[Id].image);
    ImageViewer[1]->SetImage(m_VectorData[Id].image);
    ImageViewer[2]->SetImage(m_VectorData[Id].image);
  }
  ImageViewer[0]->SetOrientation(0);
  sliceNumberSlider[0]->range( 0.0f, ImageViewer[0]->numSlices()-1 );
  ImageViewer[1]->SetOrientation(1);
  sliceNumberSlider[1]->range( 0.0f, ImageViewer[1]->numSlices()-1 );
  ImageViewer[2]->SetOrientation(2);
  sliceNumberSlider[2]->range( 0.0f, ImageViewer[2]->numSlices()-1 );
  ///////
  
  Fl::check();
  ImageViewer[0]->redraw();
  ImageViewer[1]->redraw();
  ImageViewer[2]->redraw();
  ImageViewer3D->redraw();
  Fl::check();
}

void UsimagToolConsole::ViewModeStrain3D() {
  // variable m_viewmode para saber en que modo estamos
  m_viewmode = 3;
  m_viewSplit->hide();
  m_view3D->hide();
  ImageViewer[0]->Hide();
  ImageViewer[1]->Hide();
  ImageViewer[2]->Hide();
  //ImageViewer[3]->Hide();
  ImageColorViewer->Hide();
  ImageViewerStrain3D->ResetView();
  m_viewStrain3D->resize(410,25,800,682);
  m_viewStrain3D->show();
  Fl::check();
  ImageViewerStrain3D->redraw();
  Fl::check();
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///              AQUI EMPIEZAN LOS FILTROS              /////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void UsimagToolConsole::GenericFilter( itkFilter* myfilter ) {
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	/** Is this wraper class really necessary? */
	GenericImageToImageFilter filter;
	/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
	
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &UsimagToolConsole::OnProgress );
	myfilter->AddObserver( itk::ProgressEvent(), callback );
	myfilter->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	filter.SetFilter( myfilter );
	
	filter.SetInput( m_VectorData[m_op1->value()].image );
	try{
		filter.Update();
	}
	catch(itk::ExceptionObject &excp){
		fl_alert( excp.GetDescription() );
		return;
	}
	m_VectorData.copyData( m_destino->value()-1, myfilter->GetOutput() );
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
} 

///////////// Registration /////////////////////

void UsimagToolConsole::RegisterMLDFilter(void) 
{
	if( !m_InputImage[0] || !m_InputImage[1] )
		return;
	// Create the segmentation filter
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&UsimagToolConsole::OnProgress);
	
	RegisterMLDType::Pointer regMLD = RegisterMLDType::New();
	regMLD->AddObserver( itk::ProgressEvent(), callback );
	regMLD->SetNumberOfThreads(   static_cast<unsigned int>( Threads->value() )   );
	regMLD->SetInput1( m_InputImage[0] );
	regMLD->SetInput2( m_InputImage[1] );
	regMLD->SetNLevels( static_cast<int>(reg_nlevels->value()) );
	regMLD->SetSteps( static_cast<int>(reg_steps->value()) );
	regMLD->SetUseElasticRegularization( reg_elasticreg->value() );
	regMLD->SetUseFluidRegularization( reg_fluidreg->value() );
	regMLD->SetSigmaElastic(   static_cast<int>( reg_sigmaelastic->value() )   );
	regMLD->SetSigmaFluid(   static_cast<int>( reg_sigmafluid->value() )  );
	regMLD->SetSigmaStats( reg_sigmastats->value()  );
	regMLD->SetSigmaGradient( reg_sigmagradient->value() );
	
	try{
		regMLD->Start();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData.copyData( m_destino->value()-1, regMLD->GetOutput() );
	regMLD = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


//===========================================================================
//  Funciones de control de la cámara
//===========================================================================

void UsimagToolConsole::ZoomIn()
{
//  std::cout << "m_viewmode=" << m_viewmode << std::endl;
  switch (m_viewmode)
  {
    case 0: //2D window
//      std::cout << "m_activeinput=" << m_activeinput << std::endl;
      ImageViewer[m_activeinput]->ZoomIn();
      break;
    case 1: //3D window
      ImageViewer3D->GetDefaultCamera()->Zoom(1.41);  //Value arbitrarily chosen
      ImageViewer3D->Update();
      break;
    case 3: //3D strain window
      ImageViewerStrain3D->GetDefaultCamera()->Zoom(1.41);  //Value arbitrarily chosen
      ImageViewerStrain3D->Update();
      break;
    default: //other options not yet implemented
      return;
   }

}

void UsimagToolConsole::ZoomOut()
{
//  std::cout << "m_viewmode=" << m_viewmode << std::endl;
  switch (m_viewmode)
  {
    case 0: //2D window
//      std::cout << "m_activeinput=" << m_activeinput << std::endl;
      ImageViewer[m_activeinput]->ZoomOut();
      break;
    case 1: //3D window
      ImageViewer3D->GetDefaultCamera()->Zoom(0.71);  //Value arbitrarily chosen
      ImageViewer3D->Update();
      break;
    default: //other options not yet implemented
      return;
   }

}

//===========================================================================






