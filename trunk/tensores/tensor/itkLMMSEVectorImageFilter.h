/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLMMSEVectorImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/03/29 14:53:40 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLMMSEVectorImageFilter_h
#define __itkLMMSEVectorImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkLMMSEVectorImageFilterStep.h"
#include "itkComputeRestrictedHistogram.h"
#include "itkOtsuStatistics.h"
#include "itkOtsuThreshold.h"

namespace itk
{
/** \class LMMSEVectorImageFilter
 * \brief Applies a Rician Linear Minimum Mean Square Error Filter to an Image
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 * 
 * \ingroup IntensityImageFilters
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT LMMSEVectorImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                                          InputImageType;
  typedef TOutputImage                                         OutputImageType;

  /** Standard class typedefs. */
  typedef LMMSEVectorImageFilter                               Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  typedef typename OutputImageType::Pointer                    OutputImagePointer;
  typedef typename InputImageType::Pointer                     InputImagePointer;
  typedef typename InputImageType::ConstPointer                InputImageConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( LMMSEVectorImageFilter, ImageToImageFilter );
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType                        InputPixelType;
  typedef typename OutputImageType::PixelType                       OutputPixelType;

  typedef typename InputImageType::RegionType                       InputImageRegionType;
  typedef typename OutputImageType::RegionType                      OutputImageRegionType;
  typedef typename InputImageType::SizeType                         InputSizeType;
	
  typedef LMMSEVectorImageFilterStep<TInputImage,TOutputImage>      LMMSEStepType;
  typedef typename LMMSEStepType::Pointer                           LMMSEStepPointer;

  /** Typedefs for scalar processing */
  typedef float                                                                         ScalarPixelType;
  typedef itk::Image< ScalarPixelType, TInputImage::ImageDimension >                    ScalarImageType;
  typedef itk::OtsuStatistics<InputImageType,ScalarImageType>                           LocalStatsType;
  typedef typename LocalStatsType::Pointer                                              LocalStatsPointer;
  typedef itk::OtsuThreshold<ScalarImageType, ScalarImageType>                          ThresholdType;
  typedef typename ThresholdType::Pointer                                               ThresholdPointer;
  typedef itk::ComputeRestrictedHistogram< ScalarImageType, ScalarImageType >           HistogramType;
  typedef typename HistogramType::Pointer                                               HistogramPointer;
	
  /** Set and get the radius of the neighbourhoods for estimation and filtering used to compute the mean. */
  itkSetMacro( RadiusFiltering, InputSizeType );
  itkSetMacro( RadiusEstimation, InputSizeType );
  itkGetConstReferenceMacro(RadiusFiltering, InputSizeType);
  itkGetConstReferenceMacro(RadiusEstimation, InputSizeType);

  /** Set and get the number of filter iterations. */
  itkSetMacro( Iterations, unsigned int );
  itkGetMacro( Iterations, unsigned int );

  /** Set the id for the first baseline image used to compute the noise from. */
  itkSetMacro( FirstBaseline, unsigned int );
  itkGetMacro( FirstBaseline, unsigned int );

  /** Fix behaviour in case negative values are obtained */
  itkGetMacro( UseAbsoluteValue, bool );
  itkSetMacro( UseAbsoluteValue, bool );
  itkBooleanMacro( UseAbsoluteValue );
  itkGetMacro( KeepValue, bool );
  itkSetMacro( KeepValue, bool );
  itkBooleanMacro( KeepValue );


  /** Minimum and maximum allowed noise standard deviations.*/
  itkSetMacro( MinimumNoiseSTD, double );
  itkGetMacro( MinimumNoiseSTD, double );
  itkSetMacro( MaximumNoiseSTD, double );
  itkGetMacro( MaximumNoiseSTD, double );

  /** Set and get the minimum number of voxels used for estimation and filtering. */
  itkSetMacro(MinimumNumberOfUsedVoxelsFiltering, int );
  itkGetMacro(MinimumNumberOfUsedVoxelsFiltering, int );

  itkSetMacro( Channels, unsigned int );
  itkGetConstReferenceMacro( Channels, unsigned int );

  /** This filter requires the whole input to produce its output */
  virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);	
protected:
  LMMSEVectorImageFilter();
  virtual ~LMMSEVectorImageFilter() {}
  void GenerateData();
  void PrintSelf( std::ostream& os, Indent indent) const;
private:
  LMMSEVectorImageFilter(const Self&); // purposely not implemented
  void operator=(const Self&);         // purposely not implemented

  InputSizeType      m_RadiusEstimation;
  InputSizeType      m_RadiusFiltering;
  unsigned int       m_Iterations;
  bool               m_UseAbsoluteValue;
  bool               m_KeepValue;
  int                m_MinimumNumberOfUsedVoxelsFiltering;
  double             m_MinimumNoiseSTD;
  double             m_MaximumNoiseSTD;
  unsigned int       m_FirstBaseline;
  unsigned int       m_Channels;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLMMSEVectorImageFilter.txx"
#endif

#endif
