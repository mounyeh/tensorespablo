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
// Para dar soporte de tensores: OJO!! estos includes hay que ponerlos al principio 
// porque la funcion itkImageFileReader existe en el ITK standard, y necesitamos 
// cambiarla para dar soporte de lectura de tensores
#include "tensor/itkImageFileReader.h"
#include "tensor/itkConvertPixelBuffer.h"
#include "Kretz/itkTensorVTKImageIOFactory.h"
#include "Kretz/itkTensorVTKImageIO.h"
// Usando clase de tensores de Antonio: 
#include "tensor/itkDTITensor.h"
// Usando clase de tensores de ITK:
//#include <itkDiffusionTensor3D.h>
#include "tensor/itkTractography.h"

#include <itkImage.h>
#include <itkRGBPixel.h>
#include <itkRawImageIO.h>
#include <itkMetaImageIO.h>
#include <itkJPEGImageIO.h>
#include <itkTIFFImageIO.h>
#include <itkPNGImageIO.h>
#include <itkDicomImageIO.h>
#include <itkNrrdImageIO.h>
//#include <itkNiftiImageIO.h>
#include <itkImageSeriesReader.h>
#include <itkDICOMImageIO2.h>
#include <itkDICOMSeriesFileNames.h>
//#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkVTKImageExport.h>
#include <itkExtractImageFilter.h>
#include <strain/itkComputeStrainScalars.h>

#include <itkVector.h>
#include <itkCovariantVector.h>

#include "Kretz/itkImageIOFactory.h"
#include "VtkFltk/fltkVTKImageViewer.h"

#include "Demons3D/itkMLD.h"
#include "TensorConsole.h"
#include "SegmentationConsole.h"
#include "FilteringConsole.h"
#include "BasicOpConsole.h"

#include "MyfltkImageViewer/MyfltkImageViewer.h"
#include "MyfltkImageViewer/MyfltkColorImageViewer.h"
#include "MyfltkImageViewer/GLSliceView.h"
#include "MyfltkImageViewer/GLColorSliceView.h"
#include "MyfltkImageViewer/ImageViewer.h"
#include <FL/Fl_File_Chooser.H>
#include "tensor/itkDICOMtoDWIReader.h"
#include "tensor/itkNrrdToDWIReader.h"
#include "tensor/itkDWImages.h"

#include <vector>
#include <string>
#include <iostream>
#include "VolumesContainer.h"

#include <itkNiftiImageIO.h>

#include "strain/itkStrainTensor.h"

class UsimagToolBase
{
public:

  enum { Dimension =  3 };

  //typedef unsigned char                              InputPixelType;
  typedef float                                      InputPixelType;
  typedef unsigned char                              UCharPixelType;
  typedef float                                      FloatPixelType;
  typedef itk::Image< InputPixelType, Dimension >    InputImageType;
  typedef itk::Image< InputPixelType, 2 >            InputImageType2D;
  typedef itk::Image< UCharPixelType, Dimension >    UCharImageType;
  typedef itk::Image< FloatPixelType, Dimension >    FloatImageType;
  typedef itk::Image< FloatPixelType, 2 >            FloatImageType2D;
  typedef itk::RGBPixel<float>               RGBPixelType;
  // Soporte imagenes vectoriales 
  typedef itk::Image< RGBPixelType, 3 >         ImageRGBType;
  typedef itk::Vector<float, 2> VectorPixelType;
  typedef itk::Image< VectorPixelType, 2 >    VectorImageType;
  typedef itk::CovariantVector<float, 3> Vector3DPixelType;
  typedef itk::Image< Vector3DPixelType, 3 >    Vector3DImageType;
  typedef itk::DWImages< InputPixelType, 3 >    DWImagesType;

  typedef InputImageType::RegionType                 InputRegionType;
  typedef InputImageType::RegionType::SizeType       SizeType;
  typedef InputImageType::IndexType                  InputIndexType;
  // Soporte imagenes tensoriales 
  // Usando clase de tensores de Antonio: 
  typedef itk::DTITensor<float> TensorPixelType;
  // Usando clase de tensores de ITK:
  //typedef itk::DiffusionTensor3D<float> TensorPixelType;
  typedef itk::Image<TensorPixelType, Dimension> TensorImageType;

  typedef itk::StrainTensor<float> STPixelType;
  typedef itk::Image<STPixelType, 4> STImageType;

  // Readers: 
  typedef itk::ImageFileReader< InputImageType >     FileReaderType;
  typedef itk::ImageFileReader< ImageRGBType >       ColorReaderType;
  typedef itk::ImageFileReader< FloatImageType >     FileFloatReaderType;
  typedef itk::ImageFileReader< FloatImageType2D >   FileFloat2DReaderType;
  typedef itk::ImageSeriesReader< FloatImageType >   SeriesReaderType;
  typedef itk::ImageFileReader< VectorImageType >    FileVectorReaderType;
  typedef itk::ImageFileReader< TensorImageType >    FileTensorReaderType;
  typedef itk::ImageFileReader< UCharImageType >     FileUCharReaderType;
  typedef itk::ImageFileReader< STImageType >        FileSTReaderType;
  
  // Readers IO
  typedef itk::MetaImageIO MetaImageIOType;
  typedef itk::NiftiImageIO NiftiImageIOType;
  typedef itk::JPEGImageIO JPEGImageIOType;
  typedef itk::TIFFImageIO TIFFImageIOType;
  typedef itk::PNGImageIO PNGImageIOType;
  typedef itk::TensorVTKImageIO TensorVTKImageIOType;
  typedef itk::GDCMImageIO DicomImageIOType;
  typedef itk::DICOMImageIO2 DicomSeriesImageIOType;
  typedef itk::NrrdImageIO NrrdImageIOType;
  typedef itk::RawImageIO<InputPixelType,Dimension>  RawReaderType;
  typedef itk::RawImageIO<InputPixelType,Dimension>  RawWriterType;
  typedef itk::DICOMtoDWIReader< DWImagesType >      DICOMtoDWImagesReaderType;
  typedef itk::NrrdToDWIReader< DWImagesType >       NrrdToDWImagesReaderType;

  typedef itk::RawImageIO<InputPixelType,Dimension>::IOComponentType PixelType;    
  // Writers
  typedef itk::ImageFileWriter< UCharImageType >     UCharWriterType;
  typedef itk::ImageFileWriter< FloatImageType2D >   FloatWriter2DType;
  typedef itk::ImageFileWriter< FloatImageType >     FloatWriterType;
  typedef itk::ImageFileWriter< TensorImageType >    TensorWriterType;
  typedef itk::ImageFileWriter< DWImagesType >       DWIWriterType;
  typedef itk::ImageFileWriter< STImageType >       STWriterType;

  // Viewers 
  typedef fltk::VTKImageViewer<InputImageType::PixelType>  VTKImageViewerType;
  //typedef  GLSliceView< ImagePixelType, OverlayPixelType > GLSliceViewType;
  typedef fltk::MyImageViewer<InputPixelType,UCharPixelType>     MyViewerType;
  typedef fltk::MyColorImageViewer<float,float>     MyColorViewerType;
  typedef ImageViewer<FloatPixelType,InputPixelType>             ViewerTypeFloat;
  typedef ImageViewer<InputPixelType,InputPixelType>             ViewerType;
  // Definicion de Modulos
  typedef TensorConsole TensorConsoleType;
  typedef SegmentationConsole SegmentationConsoleType;
  typedef FilteringConsole FilteringConsoleType;
  typedef BasicOpConsole BasicOpConsoleType;
  // Casters
  typedef itk::CastImageFilter<UCharImageType,FloatImageType> CastCharToRealFilterType;
  typedef itk::CastImageFilter<FloatImageType,UCharImageType> CastRealToCharFilterType;
  typedef itk::CastImageFilter<InputImageType2D,InputImageType> Cast2Dto3DFilterType;
  typedef itk::CastImageFilter<InputImageType,InputImageType2D> Cast3Dto2DFilterType;
  // Iterators
  typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
  typedef itk::ImageRegionConstIterator<ImageRGBType> ConstRGBIteratorType;
  typedef itk::ImageRegionIterator<InputImageType> IteratorType;
  typedef itk::ImageRegionIterator<ImageRGBType> RGBIteratorType;
  typedef itk::ImageRegionConstIterator<TensorImageType> ConstTensorIteratorType;
  typedef itk::ImageRegionIterator<TensorImageType> TensorIteratorType;

  //typedef DPADFilterType::TensorFilterType                     DPADDiffusionType;

  typedef itk::MLD< InputImageType, InputImageType > RegisterMLDType;
  // Connector VTK
  typedef itk::ImageToVTKImageFilter< InputImageType > ConnectorType;  
  typedef itk::VTKImageExport< InputImageType> VTKexportType;
  typedef itk::VTKImageExport< Vector3DImageType> VTKexportVectorType;
  typedef itk::VTKImageExport< InputImageType > VTKstrainExportType;
  
  typedef itk::Image< STPixelType, 3 >    				STImageType3D;
  typedef itk::ExtractImageFilter< STImageType, STImageType3D >		ExtractFilterType;
  typedef itk::ComputeStrainScalars <STImageType3D,InputImageType>     	ComputeStrainScalarsType;
  typedef ComputeStrainScalarsType::Pointer                            	ComputeStrainScalarsPointer;

	/*
  struct DataElementType {
     unsigned int Id;
	 std::string nombre;
	 InputImageType::Pointer image;
  };
  typedef std::vector<DataElementType> VectorOfDataType;
  VectorOfDataType m_VectorData;
  
  struct DataTensorElementType {
     unsigned int Id;
	 std::string nombre;
	 TensorImageType::Pointer image;
  };
  typedef std::vector<DataTensorElementType> VectorOfTensorDataType;
  VectorOfTensorDataType m_VectorTensorData;

  struct DataModelElementType {
     unsigned int Id;
	 std::string nombre;
	 vtkPolyData* data;
	 vtkActor* actor;
  };
  typedef std::vector<DataModelElementType> VectorOfModelType;
  VectorOfModelType m_VectorModelData;
	 */
	typedef VolumesContainer<DataElementType>       VectorOfDataType;
	typedef VolumesContainer<DataTensorElementType> VectorOfTensorDataType;
	typedef VolumesContainer<DataDWIElementType>    VectorOfDWIDataType;
	typedef VolumesContainer<DataModelElementType>  VectorOfModelType;
	typedef VolumesContainer<DataSTElementType>	VectorOfSTDataType;
	
	VectorOfDataType       m_VectorData;
	VectorOfTensorDataType m_VectorTensorData;
	VectorOfDWIDataType    m_VectorDWIData;
	VectorOfModelType      m_VectorModelData;
	VectorOfSTDataType     m_VectorSTData;
  
  virtual void Load(void){};
  virtual void LoadRaw(void){};
  virtual void LoadRaw(const char* filename){};
  virtual void LoadMeta(void){};
  // virtual void SaveFloat(void);
  // virtual void SaveFloat(const char * outputFileName);

  virtual void SetDimensionX( unsigned int );
  virtual void SetDimensionY( unsigned int );
  virtual void SetDimensionZ( unsigned int );
  virtual void SetSpacingX( double );
  virtual void SetSpacingY( double );
  virtual void SetSpacingZ( double );

  virtual void Show(void){};
  virtual void ShowImageSrc(int window){};
  virtual void ShowImageOut(unsigned int n, InputImageType::Pointer image);
  virtual void ShowVTK(void);

  virtual void Distribuir(InputImageType::Pointer aux, int window, const char* filename){};

  virtual void Quit(void){};
  
  unsigned int m_activeinput;
  unsigned int m_NAuxWindows;
  
  UsimagToolBase();
  virtual ~UsimagToolBase();
  
  const char * m_filename;

protected:
  
  
  // Writers 
  //FloatWriterType::Pointer    m_MetaImageWriter;
  //FloatWriterType::Pointer    m_FloatImageWriter;
  //RawReaderType::Pointer      m_RawWriter; 

  //MetaImageIOType::Pointer    m_MetaWriter;

  // Imagenes de entrada 
  InputImageType::Pointer     m_InputImage[4];
  ImageRGBType::Pointer       m_RGBImage;
  InputRegionType             m_RegionInit; 
  SizeType                    m_Size;
  FloatImageType2D::Pointer   m_FloatInputImage;

  VectorImageType::Pointer     m_VectorImage;
  TensorImageType::Pointer     m_TensorImage;

  // Imagenes de salida e intermedias
  FloatImageType::Pointer     m_FloatImage;
  FloatImageType::Pointer     m_FloatImage2;   
  InputImageType::Pointer     m_OutputImage;
  UCharImageType::Pointer     m_ImageOverlay;
  
  // Visores 
  std::vector<ViewerType>     m_VectorViewerOut;
  ViewerType                  m_v[20];
  ViewerTypeFloat             m_ViewerFloat;
  ViewerTypeFloat             m_ViewerFloat2;
  // PARA VTK 
  VTKImageViewerType::Pointer m_VTKviewer;
  ConnectorType::Pointer      m_connector;
  VTKexportType::Pointer       m_VTKexporter; 
  VTKexportVectorType::Pointer m_VTKexporterVector;
  VTKexportType::Pointer	m_VTKexporterStrain;
  
	/*
  // casters
  CastCharToRealFilterType::Pointer m_char_to_real;
  CastRealToCharFilterType::Pointer m_real_to_char;
  Cast2Dto3DFilterType::Pointer     m_2Dto3DCaster;
  Cast3Dto2DFilterType::Pointer     m_3Dto2DCaster;
  Cast2Dto3DFilterType::Pointer     m_SegmentedCaster;
  Cast3Dto2DFilterType::Pointer     m_InputCaster;
  Cast2Dto3DFilterType::Pointer     m_OutputCaster;
	 */
  
  int                         m_viewmode;
  int                         m_FileType;
  PixelType                   m_PixelType;
  unsigned int                m_ByteOrder;
  unsigned int                m_NImagesCargadas;
  unsigned int                m_NumberOfPixelsInX;
  unsigned int                m_NumberOfPixelsInY;
  unsigned int                m_NumberOfPixelsInZ;
 
  // Parametros
  unsigned int                m_NIterations;
  //MyTimeStepType              m_TimeStep;
  double                      m_Conductance;
  float m_threshold;
  double m_inmean, m_invar, m_indifmean, m_indifvar, m_inweight;
private:
  
  ViewerType                  m_ViewerSrc1;
  ViewerType                  m_ViewerSrc2;
  ViewerType                  m_ViewerSrc3;
   
  double m_SpacingX;
  double m_SpacingY;
  double m_SpacingZ;
  
};



