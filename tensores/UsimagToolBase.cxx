/*=========================================================================

  Program:   UsimagToolBase
  Language:  C++
  Date:      5-07-2007
  Version:   1.0

  Copyright (c) 2007 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/

#include <fstream>
#include "UsimagToolBase.h"
#include "FL/fl_ask.H"

UsimagToolBase
::UsimagToolBase()
{

  //Readers
  //m_JPEGFloatImageReader = FileFloat2DReaderType::New();

  //Writers 
  //m_MetaImageWriter   =  FloatWriterType::New(); 
  //m_FloatImageWriter  =  FloatWriterType::New(); 
  //m_RawWriter         =  RawReaderType::New(); 
  //m_MetaWriter        =  MetaImageIOType::New(); 
  
  //m_MetaImageWriter->SetImageIO( m_MetaWriter );
  
  /// PARA VTK
  m_VTKviewer         = VTKImageViewerType::New();
  m_connector         = ConnectorType::New();
  //connector->GetInput(VTKreader->GetOutput);

  m_VTKexporter = VTKexportType::New(); 
  m_VTKexporterVector = VTKexportVectorType::New(); 
  
  //Parametros generales
  m_NImagesCargadas = 0;
  m_activeinput = 0;
  m_NAuxWindows = 0;
  m_viewmode    = 0;
  m_NumberOfPixelsInX = 256;
  m_NumberOfPixelsInY = 256;
  m_NumberOfPixelsInZ = 1;

  m_SpacingX = 1;
  m_SpacingY = 1;
  m_SpacingZ = 1;

  //Parametros especificos
  m_NIterations = 10;
  //m_TimeStep = 1;
  m_Conductance = 0.1;
  // Segmentacion 
  m_threshold = 100;
  // Imagenes
  m_InputImage[0] = InputImageType::New();
  m_InputImage[1] = InputImageType::New();
  m_InputImage[2] = InputImageType::New();
  m_InputImage[3] = InputImageType::New();
  m_RGBImage = ImageRGBType::New();

  m_FloatInputImage = FloatImageType2D::New();
  m_VectorImage = VectorImageType::New();
  m_TensorImage = TensorImageType::New();
  m_ImageOverlay = UCharImageType::New();
	
  m_Size[0] = 128; m_Size[1] = 128; m_Size[2] = 1; 
  m_RegionInit.SetSize(m_Size);
  for (int i=0; i < 4; i++ ) {
    m_InputImage[i]->SetRegions(m_RegionInit);
    m_InputImage[i]->Allocate();
    m_InputImage[i]->FillBuffer(0);
  }
  m_RGBImage->SetRegions(m_RegionInit);
  m_RGBImage->Allocate();
  RGBPixelType RGBpixel;
  RGBpixel[0] = 0; RGBpixel[1] = 0; RGBpixel[2] = 0; 
  m_RGBImage->FillBuffer(RGBpixel);

  m_OutputImage = InputImageType::New();
  m_FloatImage  = FloatImageType::New();
  m_FloatImage2  = FloatImageType::New();

  //// Para reservar memoria para una imagen ////
  //ImageType::SizeType imSize = {{256, 256, 10}};
  //ImageType::Pointer im = ImageType::New();
  //im->SetRegions(imSize);
  //im->Allocate();

}

UsimagToolBase::~UsimagToolBase()
{
}

void UsimagToolBase::SetDimensionX( unsigned int numberOfPixels )
{
  m_NumberOfPixelsInX = numberOfPixels;
}


void UsimagToolBase::SetDimensionY( unsigned int numberOfPixels )
{
  m_NumberOfPixelsInY = numberOfPixels;
}


void UsimagToolBase::SetDimensionZ( unsigned int numberOfPixels )
{
  m_NumberOfPixelsInZ = numberOfPixels;
}

void UsimagToolBase::SetSpacingX( double value )
{
  m_SpacingX = value;
}

void UsimagToolBase::SetSpacingY( double value )
{
  m_SpacingY = value;
}

void UsimagToolBase::SetSpacingZ( double value )
{
  m_SpacingZ = value;
}


void UsimagToolBase::ShowImageOut(unsigned int n, InputImageType::Pointer image)
{
  if (n < 10) {
    char label[50];
    sprintf(label,"Aux %d",n);
    m_v[n].SetLabel(label);
    m_v[n].SetImage( image );
    m_v[n].Show();
	m_NAuxWindows++;
  }	
}


void UsimagToolBase::ShowVTK(void)
{
  m_VTKviewer->SetImage(m_InputImage[0]);
  m_VTKviewer->Show();
}



/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
/** TO DO: Are these methods necessary? */

/*
void UsimagToolBase::SaveFloat( void )
{
  
  m_filename = fl_file_chooser("Image Filename","*","");
 
  if( !m_filename )
    {
    return;
    }

  this->SaveFloat( m_filename );
}
*/

/*
void UsimagToolBase::SaveFloat( const char * filename )
{
	// m_FloatImageWriter->SetFileName( filename );
	// m_FloatImageWriter->SetInput( m_DicomSeriesImageReader->GetOutput());
	
	// Attempt to Write
	try{
		// m_FloatImageWriter->Write();
	}
	catch (itk::ExceptionObject &ex){
		std::cout << ex << std::endl;
	}
}
*/

/** $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */


