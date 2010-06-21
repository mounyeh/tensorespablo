/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkOtsuStatistics.h,v $
  Language:  C++
  Date:      $Date: 2008/02/7 14:28:51 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkThresholdToMaskImageFilter_h
#define __itkThresholdToMaskImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{
/** \class OtsuThreshold */
	
template <class TInputImage, class TOutputImage>
class ITK_EXPORT ThresholdToMaskImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                               InputImageType;
  typedef typename InputImageType::Pointer          InputImagePointer;
  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef TOutputImage                              OutputImageType;
  typedef typename OutputImageType::Pointer         OutputImagePointer;
  typedef typename OutputImageType::ConstPointer    OutputImageConstPointer;
	
  /** Standard class typedefs. */
  typedef ThresholdToMaskImageFilter                           Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ThresholdToMaskImageFilter, ImageToImageFilter );
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType    InputSizeType;

  /** Set and get the number of DWI channels. */
  itkGetMacro( Threshold, double );
  itkSetMacro( Threshold, double );
protected:
  ThresholdToMaskImageFilter();
  virtual ~ThresholdToMaskImageFilter() {}
  // Threaded filter!
  void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
  ThresholdToMaskImageFilter(const Self&);   // purposely not implemented
  void operator=(const Self&);               // purposely not implemented
  // Final Threshold
  double m_Threshold;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkThresholdToMaskImageFilter.txx"
#endif

#endif
