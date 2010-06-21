/*=========================================================================

  Program:   BasicOpConsole.h
  Language:  C++
  Date:      7-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#ifndef _BasicOpConsole_h
#define _BasicOpConsole_h

#include "BasicOpGUI.h"
#include "Viewer3D.h"
// fltk includes
#include <FL/Fl_Choice.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Value_Output.H>
#include <fltkProgressBar.h>
// itk includes
#include <itkChangeLabelImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkErodeObjectMorphologyImageFilter.h>
#include <itkDilateObjectMorphologyImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryCrossStructuringElement.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkInvertIntensityImageFilter.h>
#include <itkRGBPixel.h>
// Other includes
#include "MyfltkImageViewer/MyfltkImageViewer.h"
#include "MyfltkImageViewer/MyfltkColorImageViewer.h"
#include "VolumesContainer.h"


class BasicOpConsole : public BasicOpGUI
{
public:
  enum { Dimension =  3 };
  typedef float                                      InputPixelType;
  typedef float                                      FloatPixelType;
  typedef itk::Image< FloatPixelType, Dimension >    FloatImageType;
  // Usando clase de tensores de ITK:
  typedef itk::ImageRegionIterator<InputImageType> IteratorType;
  typedef itk::RGBPixel<float>             RGBPixelType;
  typedef itk::Image< RGBPixelType, 3 >            ImageRGBType;
  typedef itk::ImageRegionIterator<ImageRGBType>   RGBIteratorType;
  typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
  typedef fltk::MyColorImageViewer<float,float>     MyColorViewerType;
  
  // Filtros:
  typedef itk::GradientMagnitudeImageFilter< InputImageType,InputImageType  >  GradientFilterType;
  typedef itk::ChangeLabelImageFilter<InputImageType,InputImageType> RelabelFilterType;
  typedef itk::RescaleIntensityImageFilter<InputImageType,InputImageType> RescaleFilterType;
  typedef itk::RelabelComponentImageFilter<InputImageType,InputImageType> RelabelCompFilterType;
  typedef itk::AddImageFilter<InputImageType,InputImageType,InputImageType> AddFilterType;
  typedef itk::MultiplyImageFilter<InputImageType,InputImageType,InputImageType> MultiplyFilterType;
  typedef itk::InvertIntensityImageFilter< InputImageType, InputImageType> InvertFilterType;
  // Structured Element 
  typedef itk::BinaryBallStructuringElement< InputPixelType, Dimension> StructuringElementType;
  //typedef itk::BinaryCrossStructuringElement< InputPixelType, Dimension> StructuringElementType;
  typedef itk::ErodeObjectMorphologyImageFilter<InputImageType, InputImageType, StructuringElementType > ErodeFilterType;
  typedef itk::DilateObjectMorphologyImageFilter<InputImageType, InputImageType, StructuringElementType > DilateFilterType;
  // Definition of types for progress bar
  typedef itk::MemberCommand<BasicOpConsole> CommandType;
  typedef itk::SmartPointer<CommandType> CommandPointer;
  typedef fltk::ProgressBar ProgressBarType;

  typedef VolumesContainer<DataModelElementType> VectorOfModelType;
  typedef VolumesContainer<DataElementType>      VectorOfDataType;

  // Funciones miembro de la clase
  BasicOpConsole(int X, int Y, int W, int H, const char *L = 0);
  virtual ~BasicOpConsole();

  //void SetActiveInput(unsigned int &activeinput);
  void SetModelDataBrowser(Fl_Browser *modeldataBrowserIn);
  void SetImageViewer3D(Viewer3D* ImageViewer3DIn);
  void SetImageColorViewer(MyColorViewerType* ImageColorViewer3DIn);
  void SetVectorModelData(void* VectorModelData);
  void SetVectorData(void* VectorData);
  void SetSliders(Fl_Value_Slider **sliders); 
  void SetThreads(int threads); 
  void SetProgressSliders(ProgressBarType *progressbar); 
  void SetProgressCounter(Fl_Value_Output *progressCounter);
  void OnProgress(itk::Object *object, const itk::EventObject &event);
	
  void AddFilter(void);
  void MultiplyFilter(void);
  void RelabelFilter(void);
  void RescaleFilter(void);
  void ErodeFilter(void);
  void DilateFilter(void);
  void RelabelCompFilter(void);
  void InvertFilter(void);
  void GradientFilter(void);
  
  VectorOfModelType*           m_VectorModelData;
  VectorOfDataType*            m_VectorData;
  Fl_Value_Output*             m_ProgressCounter;
  ProgressBarType*             progressSlider;
  MyColorViewerType*           ImageColorViewer;
  Viewer3D*                    ImageViewer3D;
  Fl_Browser*                  m_modeldataBrowser;
  Fl_Value_Slider**            sliceNumberSlider;
  int                          m_threads;
  int                          m_NImagesCargadas;
};

#endif




