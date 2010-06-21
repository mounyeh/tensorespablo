/*=========================================================================

  Program:   SegmentationConsole.h
  Language:  C++
  Date:      27-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#ifndef _SegmentationConsole_h
#define _SegmentationConsole_h

#include "SegmentationGUI.h"

//fltk includes:
#include <FL/Fl_Choice.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Value_Output.H>
#include <fltkProgressBar.h>
// Viewer inckudes
#include "MyfltkImageViewer/MyfltkImageViewer.h"
#include "MyfltkImageViewer/MyfltkColorImageViewer.h"
#include "Viewer3D.h"
// itk includes:
#include <itkRGBPixel.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkThresholdSegmentationLevelSetImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkFastMarchingImageFilter.h>
#include <itkWatershedImageFilter.h>
#include <itkVTKImageExport.h>
#include <itkCastImageFilter.h>
// My fiter includes: 
#include "knn-1canal/classifyKNNCore.h"
#include "knn-1canal/voronoiFilter.h"
#include "vtkMarcacionElipse/vtkMarcacionElipsoide.h"
#include "vtkMarcacionElipse/vtkPlantillaAjustada.h"
#include "vtkMarcacionElipse/vtkFuncionVerosimilitud.h"
#include "vtkMarcacionElipse/vtkOptimizaContorno.h"
#include "vtkMarcacionElipse/vtkValidacion.h"
#include "VolumesContainer.h"

namespace itk {
  template<class TObject> class MemberCommand;
  class EventObject;
}

class SegmentationConsole : public SegmentationGUI
{
public:
  enum { Dimension =  3 };
  // Definition of pixel types and images
  typedef float                                      InputPixelType;
  typedef float                                      FloatPixelType;
  typedef itk::Image< FloatPixelType, Dimension >    FloatImageType;
  typedef itk::Image< FloatPixelType, 2 >            FloatImageType2D;
  typedef itk::Image< InputPixelType, 2 >            InputImageType2D;
  typedef itk::RGBPixel<float>                     RGBPixelType;
  typedef itk::Image< RGBPixelType, 3 >            ImageRGBType;
  //typedef itk::Image< InputPixelType, Dimension >    InputImageType;
  // Definitions of types for KNN segmentation:
  typedef itk::Image<int, 1> VoronoiInputType;
  typedef VoronoiInputType::IndexType IndexType;
  typedef itk::voronoiFilter<VoronoiInputType, VoronoiInputType::PixelType, 1> VoronoiFilterType;
  typedef itk::Image< std::vector<int> , 1>  DiagramImageType; 
  typedef itk::classifyKNNCore<InputImageType, InputImageType, 3> classifyKNNCoreType;  
  // Definitions of types for segmentation filters
  typedef itk::ThresholdSegmentationLevelSetImageFilter<FloatImageType2D,FloatImageType2D> SegmentationFilterType;
  typedef itk::BinaryThresholdImageFilter<FloatImageType2D,InputImageType2D >   BinaryThresholdType;
  typedef itk::BinaryThresholdImageFilter<InputImageType, InputImageType> BinaryThresholdFilterType;
  typedef itk::OtsuThresholdImageFilter<InputImageType, InputImageType> OtsuThresholdFilterType;
  typedef itk::WatershedImageFilter<InputImageType> WatershedFilterType;
  typedef itk::FastMarchingImageFilter<FloatImageType2D, FloatImageType2D >   FastMarchingFilter2DType;
  typedef itk::FastMarchingImageFilter<FloatImageType, FloatImageType >   FastMarchingFilter3DType;
  typedef itk::VTKImageExport< InputImageType> VTKexportType;
  typedef itk::CastImageFilter<InputImageType2D,InputImageType> Cast2Dto3DFilterType;
  typedef itk::CastImageFilter<InputImageType,InputImageType2D> Cast3Dto2DFilterType;
  typedef FastMarchingFilter2DType::NodeType       NodeType2D;
  typedef FastMarchingFilter2DType::NodeContainer  NodeContainer2D;
  // Definition of iterators and viewers
  typedef itk::ImageRegionIterator<InputImageType>      IteratorType;
  typedef itk::ImageRegionIterator<ImageRGBType>        RGBIteratorType;
  typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
  typedef fltk::MyColorImageViewer<float,float>         MyColorViewerType;
  // Definition of types for progress bar
  typedef itk::MemberCommand<SegmentationConsole> CommandType;
  typedef itk::SmartPointer<CommandType> CommandPointer;
  typedef fltk::ProgressBar ProgressBarType;
  
  
  // Member functions
  SegmentationConsole(int X, int Y, int W, int H, const char *L = 0);
  virtual ~SegmentationConsole();

  void SetModelDataBrowser(Fl_Browser *modeldataBrowserIn);
  void SetImageViewer3D(Viewer3D* ImageViewer3DIn);
  void SetImageColorViewer(MyColorViewerType* ImageColorViewer3DIn);
  void SetSliders(Fl_Value_Slider **sliders);
  void SetProgressSliders(ProgressBarType *progressbar); 
  void SetProgressCounter(Fl_Value_Output *progressCounter);
  void SetThreads(int threads); 
  void SetNImagesCargadas(unsigned int &NImagesCargadas); 
  void OnProgress(itk::Object *object, const itk::EventObject &event);
  void GenerateNewData(InputImageType::Pointer aux);
  void CopyData(InputImageType::Pointer aux);
  
  void GetClickPoints(float x, float y);
  static void ClickSelectCallback(float x, float y, float z, float v, int c, void * args );

  void LoadME();
  void LoadKNN();
  void KNNFileInputChange(void);
  void KnnFilter();
  
  void MarcacionElipseInit(void);
  void CalculoModeloVerosimilitud(void);
  void CalculoOptimizacionSA(void);

  void ThresholdFilter(void);
  void OtsuThresholdFilter(void);
  void ResetAllParameters(void);
  void UpdateViewerAfterIteration(void);
  void ThresholdLevelSetSegmentation(void);
  
	typedef VolumesContainer<DataModelElementType> VectorOfModelType;
	typedef VolumesContainer<DataElementType>      VectorOfDataType;
	VectorOfModelType*                             m_VectorModelData;
	VectorOfDataType*                              m_VectorData;
	void SetVectorModelData(void* VectorModelData);
	void SetVectorData(void* VectorData);
	
  // Member atributes
  vtkMarcacionElipsoide* m_vtkMarcacionElipsoide;
  vtkPlantillaAjustada* m_vtkPlantillaAjustada;
  vtkFuncionVerosimilitud* m_vtkFuncionVerosimilitud;
  vtkOptimizaContorno* m_vtkOptimizaContorno;
  vtkValidacion* m_vtkValidacion;
  
  SegmentationFilterType::Pointer    m_LevelSetSegmentation; 
  BinaryThresholdType::Pointer       m_thresholder;    
  WatershedFilterType::Pointer       m_WatershedFilter;
  VTKexportType::Pointer             m_VTKexporter; 
  Cast2Dto3DFilterType::Pointer      m_SegmentedCaster;
  Cast3Dto2DFilterType::Pointer      m_InputCaster;
  Cast2Dto3DFilterType::Pointer      m_OutputCaster;
  FloatImageType2D::Pointer          m_FloatInputImage;
  
  FastMarchingFilter2DType::Pointer  m_fastMarching;
  FloatImageType2D::IndexType        m_seedPosition2D;
  NodeType2D                         m_node2D;
  NodeContainer2D::Pointer           m_seeds2D;
  
  
  MyColorViewerType* ImageColorViewer;
  Viewer3D* ImageViewer3D;
  Fl_Browser *m_modeldataBrowser;
  Fl_Value_Slider **sliceNumberSlider;
  Fl_Value_Output *m_ProgressCounter;
  ProgressBarType *progressSlider;
  int m_threads,m_NImagesCargadas;
  const char * ME_filename;
  const char * KNN_filename;
   
};

#endif









