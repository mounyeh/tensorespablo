/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLocalMeanFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLocalMeanFilter_h
#define __itkLocalMeanFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

#include "itkRecursiveGaussianImageFilter.h"

namespace itk
{
/** \class LocalMeanFilter
 * 
 */
template < class TScalar, unsigned int NDimension >
class ITK_EXPORT LocalMeanFilter :
	public ImageToImageFilter<    itk::Image< TScalar, NDimension >,   itk::Image< TScalar, NDimension >    >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension,  unsigned int, NDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef itk::Image< TScalar, NDimension >                                 InputImageType;
  typedef itk::Image< TScalar, NDimension >                                 OutputImageType;

  /** Standard class typedefs. */
  typedef LocalMeanFilter                                                   Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >             Superclass;
  typedef SmartPointer< Self >                                              Pointer;
  typedef SmartPointer< const Self >                                        ConstPointer;
  
  typedef itk::RecursiveGaussianImageFilter< InputImageType, OutputImageType >
		                                                                    SmoothingFilterType;
  typedef typename SmoothingFilterType::Pointer                             SmoothingFilterPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(  Self  );

  /** Run-time type information (and related methods). */
  itkTypeMacro( LocalMeanFilter, ImageToImageFilter );
  
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
  itkGetMacro( Sigma, double );
  itkSetMacro( Sigma, double );

protected:
  LocalMeanFilter();
  virtual ~LocalMeanFilter() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;
  /** Multi-threaded filter */
  void GenerateData( void );

private:
  LocalMeanFilter(const Self&);   // purposely not implemented
  void operator=(const Self&);    // purposely not implemented
  double m_Sigma;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLocalMeanFilter.txx"
#endif

#endif
