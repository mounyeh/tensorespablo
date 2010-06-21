/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeCoefficientOfVariation.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComputeCoefficientOfVariation_h
#define __itkComputeCoefficientOfVariation_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ComputeCoefficientOfVariation : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                        InputImageType;
  typedef TOutputImage                       OutputImageType;
  typedef typename TInputImage::Pointer      InputImagePointer;
  typedef typename TInputImage::ConstPointer InputImageConstPointer;
  typedef typename TOutputImage::Pointer     OutputImagePointer;

  /** Standard class typedefs. */
  typedef ComputeCoefficientOfVariation                        Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ComputeCoefficientOfVariation, ImageToImageFilter);
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType    InputSizeType;
  
  /** Set the first input */
  void SetInput1( const InputImageType  * image ){ this->SetInput( 0, image ); }
  /** Set the second input */
  void SetInput2( const InputImageType  * image ){ this->SetInput( 1, image ); }
  /** Get the first input */
  const InputImageType  * GetInput1( void ){ return this->GetInput( 0 ); }
  /** Get the second input */
  const InputImageType  * GetInput2( void ){ return this->GetInput( 1 ); }

  itkSetMacro(               NumberOfPixels, unsigned int );
  itkGetConstReferenceMacro( NumberOfPixels, unsigned int );
protected:
  ComputeCoefficientOfVariation();
  virtual ~ComputeCoefficientOfVariation() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
  void GenerateInputRequestedRegion();
private:
  ComputeCoefficientOfVariation(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  unsigned int m_NumberOfPixels;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeCoefficientOfVariation.txx"
#endif

#endif
