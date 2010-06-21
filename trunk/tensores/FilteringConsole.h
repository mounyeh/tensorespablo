/*=========================================================================

  Program:   FilteringConsole.h
  Language:  C++
  Date:      27-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#ifndef _FilteringConsole_h
#define _FilteringConsole_h

#include "FilteringGUI.h"

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
//#include <itkAnisotropicDiffusionImageFilter.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkBilateralImageFilter.h>
#include <itkAbsImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkMeanImageFilter.h>
#include <itkZeroCrossingBasedEdgeDetectionImageFilter.h>
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkVTKImageExport.h>
#include <itkCastImageFilter.h>
// My fiter includes: 
#include "wiener/itkWienerFilter.h"
#include "ASR/itkASRFilter.h"
#include "SRAD/itkSRADFilter.h"
#include "DPAD/itkDPADFilter.h"
#include "VolumesContainer.h"

namespace itk {
  template<class TObject> class MemberCommand;
  class EventObject;
}

class FilteringConsole : public FilteringGUI
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
  // Definitions of types for segmentation filters
  typedef itk::AnisotropicDiffusionImageFilter<InputImageType,InputImageType>::TimeStepType MyTimeStepType;
  typedef itk::BilateralImageFilter<InputImageType,InputImageType  >  BilateralFilterType;
  typedef itk::GradientAnisotropicDiffusionImageFilter<InputImageType,InputImageType  >  GradientDiffusionType;
  typedef itk::CurvatureAnisotropicDiffusionImageFilter<FloatImageType,FloatImageType  >  CurvatureDiffusionType;
  typedef itk::CurvatureAnisotropicDiffusionImageFilter<InputImageType2D,InputImageType2D  >  Curvature2DDiffusionType;
  typedef itk::RecursiveGaussianImageFilter<InputImageType, InputImageType > GaussianFilterType;
  typedef itk::MedianImageFilter<InputImageType, InputImageType > MedianFilterType;
  typedef itk::MeanImageFilter<InputImageType, InputImageType > MeanFilterType;
  typedef itk::AbsoluteValueDifferenceImageFilter< InputImageType, InputImageType, InputImageType> ABSValDifFilterType;
  typedef itk::ZeroCrossingBasedEdgeDetectionImageFilter<InputImageType,InputImageType> ZeroEdgeFilterType;
  typedef itk::CannyEdgeDetectionImageFilter< InputImageType, InputImageType > CannyFilterType;
  typedef itk::AbsImageFilter< InputImageType, InputImageType> ABSFilterType;
  // Casters
  typedef itk::CastImageFilter<InputImageType2D,InputImageType> Cast2Dto3DFilterType;
  typedef itk::CastImageFilter<InputImageType,InputImageType2D> Cast3Dto2DFilterType;
  // Filtros propios
  typedef itk::WienerFilter< InputImageType, InputImageType >  WienerFilterType;
  typedef itk::ASRFilter< InputImageType, InputImageType >     ASRFilterType;
  typedef itk::ASRFilter< InputImageType2D, InputImageType2D > ASRFilterType2D;
  typedef itk::SRADFilter< InputImageType, InputImageType >    SRADFilterType;
  typedef SRADFilterType::TensorFilterType                     SRADDiffusionType;
  typedef itk::DPADFilter< InputImageType, InputImageType >    DPADFilterType;
  // Definition of iterators and viewers
  typedef itk::ImageRegionIterator<InputImageType>      IteratorType;
  typedef itk::ImageRegionIterator<ImageRGBType>        RGBIteratorType;
  typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
  typedef fltk::MyColorImageViewer<float,float>         MyColorViewerType;
  // Definition of types for progress bar
  typedef itk::MemberCommand<FilteringConsole> CommandType;
  typedef itk::SmartPointer<CommandType> CommandPointer;
  typedef fltk::ProgressBar ProgressBarType;
  
  // Member functions
  FilteringConsole(int X, int Y, int W, int H, const char *L = 0);
  virtual ~FilteringConsole();

  void SetModelDataBrowser(Fl_Browser *modeldataBrowserIn);
  void SetImageViewer3D(Viewer3D* ImageViewer3DIn);
  void SetImageColorViewer(MyColorViewerType* ImageColorViewer3DIn);
  void SetSliders(Fl_Value_Slider **sliders);
  void SetProgressSliders(ProgressBarType *progressbar); 
  void SetProgressCounter(Fl_Value_Output *progressCounter);
  void SetThreads(int threads); 
  void SetNImagesCargadas(unsigned int &NImagesCargadas); 
  void OnProgress(itk::Object *object, const itk::EventObject &event);
  void OnTipoFiltradoChange();
  void OnTipoFiltroWienerChange();
  
  void ABSFilter(void);
  void ABSValDifFilter(void);
  void BilateralFilter(void);
  void CurvatureAnisotropicFilter(void); 
  void GradientDiffusionFilter(void); 
  void GaussianFilter(void);
  void MedianFilter(void);
  void MeanFilter(void);
  void WienerFilter(void);
  void ASRFilter(void);
  void SRADFilter(void);
  void DPADFilter(void);
  void CannyFilter(void);
  void ZeroEdgeFilter(void);
	
	typedef VolumesContainer<DataModelElementType> VectorOfModelType;
	typedef VolumesContainer<DataElementType>      VectorOfDataType;
	VectorOfModelType* m_VectorModelData;
	VectorOfDataType*  m_VectorData;
	void SetVectorModelData(void* VectorModelData);
	void SetVectorData(void* VectorData);
	
  MyColorViewerType* ImageColorViewer;
  Viewer3D* ImageViewer3D;
  Fl_Browser *m_modeldataBrowser;
  Fl_Value_Slider **sliceNumberSlider;
  Fl_Value_Output *m_ProgressCounter;
  ProgressBarType *progressSlider;
  int m_threads,m_NImagesCargadas;
  int m_filtroWiener;
  int m_tipofiltrado;
};

#endif




