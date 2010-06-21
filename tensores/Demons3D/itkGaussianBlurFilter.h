/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianBlurFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianBlurFilter_h
#define __itkGaussianBlurFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{
/** \class GaussianBlurFilter
 * 
 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT GaussianBlurFilter : public ImageToImageFilter<  TInputImage, TOutputImage  >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension,  unsigned int, TInputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                                                       InputImageType;
  typedef TOutputImage                                                      OutputImageType;

  /** Standard class typedefs. */
  typedef GaussianBlurFilter                                                Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >             Superclass;
  typedef SmartPointer<Self>                                                Pointer;
  typedef SmartPointer<const Self>                                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(  Self  );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GaussianBlurFilter, ImageToImageFilter );
  
  /** Image typedef support. */
  typedef typename InputImageType::Pointer          InputImagePointer;
  typedef typename OutputImageType::Pointer         OutputImagePointer;

  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef typename OutputImageType::ConstPointer    OutputImageConstPointer;

  typedef typename InputImageType::PixelType        InputPixelType;
  typedef typename OutputImageType::PixelType       OutputPixelType;

  typedef typename InputImageType::RegionType       InputRegionType;
  typedef typename OutputImageType::RegionType      OutputRegionType;

  typedef typename InputImageType::SizeType         SizeType;
  typedef typename InputImageType::IndexType        IndexType;
  typedef typename InputImageType::SpacingType      SpacingType;

  typedef typename InputPixelType::ValueType        ValueType;

  /** Set and get the size of the gaussian operator: */
  itkSetMacro(               Radius,    unsigned int );
  itkGetConstReferenceMacro( Radius,    unsigned int );
  /** Set and get the direction of the operator: */
  itkSetMacro(               Direction, unsigned int );
  itkGetConstReferenceMacro( Direction, unsigned int );
  /** Set and get the tolerance of the outlier rejection: */
  itkSetMacro(               Tol,       double       );
  itkGetConstReferenceMacro( Tol,       double       );

protected:
  GaussianBlurFilter();
  virtual ~GaussianBlurFilter() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;
  /** Multi-threaded filter */
  void ThreadedGenerateData( const OutputRegionType& outputRegionForThread, int threadId );
  /** Pre-compute the gaussian kernel */
  void BeforeThreadedGenerateData( void );
  /** We need a larger input than the output */
  virtual void GenerateInputRequestedRegion( );

private:
  GaussianBlurFilter(const Self&); //purposely not implemented
  void operator=(const Self&);     //purposely not implemented
  
  unsigned int            m_Radius;
  unsigned int            m_Direction;
  ValueType               m_Tol;
  itk::Array< ValueType > m_Operator;
  ValueType               m_Correction;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianBlurFilter.txx"
#endif

#endif
