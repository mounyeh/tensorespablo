/*=========================================================================

  Program:   UsimagToolConsoleSegmentation
  Language:  C++
  Date:      5-07-2007
  Version:   1.0

  Copyright (c) 2007 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#include <itkObject.h>
#include <string>
#include <itkGradientMagnitudeImageFilter.h>

class UsimagToolConsoleSegmentation {
public:
  enum { Dimension =  3 };
  typedef float                                      InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >    InputImageType;
  typedef itk::GradientMagnitudeImageFilter< InputImageType,InputImageType  >  GradientFilterType;

  UsimagToolConsoleSegmentation();
  virtual ~UsimagToolConsoleSegmentation();
  
  void GradientFilter(unsigned int, InputImageType*);
  void MEFileInputChange(void);
  void MarcacionElipseInit(void);
  void CalculoModeloVerosimilitud(void);
  void CalculoOptimizacionSA(void);

protected:

 GradientFilterType::Pointer        m_GradientFilter;

};

