/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianSmoothFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianSmoothFilter_h
#define __itkGaussianSmoothFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkGaussianBlurFilter.h"

namespace itk
{
/** \class GaussianSmoothFilter
 * 
 */
template < class TScalar, unsigned int NDimension >
class ITK_EXPORT GaussianSmoothFilter :
	public ImageToImageFilter<  itk::Image< itk::Point<TScalar,NDimension>, NDimension >, itk::Image< itk::Point<TScalar,NDimension>, NDimension >  >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension,  unsigned int, NDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef itk::Image< itk::Point<TScalar,NDimension>, NDimension >          InputImageType;
  typedef itk::Image< itk::Point<TScalar,NDimension>, NDimension >          OutputImageType;

  /** Standard class typedefs. */
  typedef GaussianSmoothFilter                                              Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >             Superclass;
  typedef SmartPointer<Self>                                                Pointer;
  typedef SmartPointer<const Self>                                          ConstPointer;

  typedef itk::GaussianBlurFilter< InputImageType, OutputImageType >        SmoothingFilterType;
  typedef typename SmoothingFilterType::Pointer                             SmoothingFilterPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(  Self  );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GaussianSmoothFilter, ImageToImageFilter );
  
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

  /** Set/Get the variance of the smoothing kernel */   
  itkGetMacro( Sigma, SizeType );
  itkSetMacro( Sigma, SizeType );

protected:
  GaussianSmoothFilter();
  virtual ~GaussianSmoothFilter() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;
  /** Multi-threaded filter */
  void GenerateData( void );

private:
  GaussianSmoothFilter(const Self&);   //purposely not implemented
  void operator=(const Self&);         //purposely not implemented
  SizeType m_Sigma;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianSmoothFilter.txx"
#endif

#endif
