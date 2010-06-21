/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCoefficientOfVariationToDiffusionTerm.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCoefficientOfVariationToDiffusionTerm_h
#define __itkCoefficientOfVariationToDiffusionTerm_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

/** Typedef to check the way we compute the diffusion term. */
typedef enum DifussionTermMode{__Aja__,__Yu__,__Simplified__} DiffusionTermMode;

template <class TInputImage, class TOutputImage>
class ITK_EXPORT CoefficientOfVariationToDiffusionTerm : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                        InputImageType;
  typedef TOutputImage                       OutputImageType;
  typedef typename TInputImage::ConstPointer InputImageConstPointer;
  typedef typename TOutputImage::Pointer     OutputImagePointer;

  /** Standard class typedefs. */
  typedef CoefficientOfVariationToDiffusionTerm                Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CoefficientOfVariationToDiffusionTerm, ImageToImageFilter);
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType    InputSizeType;
  
  void SetUseModeDiffusionTerm( DiffusionTermMode mode )
  {
  	m_DiffusionTermMode = mode;
  }
  DiffusionTermMode GetUseModeDiffusionTerm( void ){
  	return m_DiffusionTermMode;
  }
  void SetUseAjaDiffusionTerm( void )
  {
  	m_DiffusionTermMode = __Aja__;
  }
  void SetUseYuDiffusionTerm( void )
  {
  	m_DiffusionTermMode = __Yu__;
  }
  void SetUseSimplifiedDiffusionTerm( void )
  {
  	m_DiffusionTermMode = __Simplified__;
  }

  /** Set the radius of the neighborhood used to compute the NDSticks. */
  itkSetMacro(               Noise, double );
  itkGetConstReferenceMacro( Noise, double );
protected:
  CoefficientOfVariationToDiffusionTerm();
  virtual ~CoefficientOfVariationToDiffusionTerm() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
  CoefficientOfVariationToDiffusionTerm(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  double            m_Noise;
  DiffusionTermMode m_DiffusionTermMode;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCoefficientOfVariationToDiffusionTerm.txx"
#endif

#endif
