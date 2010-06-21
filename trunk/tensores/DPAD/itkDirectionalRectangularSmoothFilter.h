/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDirectionalRectangularSmoothFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDirectionalRectangularSmoothFilter_h
#define __itkDirectionalRectangularSmoothFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT DirectionalRectangularSmoothFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef DirectionalRectangularSmoothFilter             Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DirectionalRectangularSmoothFilter, ImageToImageFilter);

  /** Convenient typedefs for simplifying declarations. */
  typedef          TInputImage                   InputImageType;
  typedef typename InputImageType::ConstPointer  InputImageConstPointer;
  typedef typename InputImageType::Pointer       InputImagePointer;
  typedef          TOutputImage                  OutputImageType;
  typedef typename OutputImageType::Pointer      OutputImagePointer;

  /** Image typedef support. */
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  
  typedef typename InputImageType::SizeType    InputSizeType;
  typedef typename InputImageType::IndexType   InputIndexType;

  /** Set and get the radius: */
  itkSetMacro(               Radius,    unsigned int );
  itkGetConstReferenceMacro( Radius,    unsigned int );

  /** Set and get the direcion: */
  itkSetMacro(               Direction, unsigned int );
  itkGetConstReferenceMacro( Direction, unsigned int );
 
  virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);

protected:
  DirectionalRectangularSmoothFilter();
  virtual ~DirectionalRectangularSmoothFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void ThreadedGenerateData(  const OutputImageRegionType& outputRegionForThread, int threadId );

private:
  DirectionalRectangularSmoothFilter(const Self&); //purposely not implemented
  void operator=(const Self&);                     //purposely not implemented
  // Parameters:
  unsigned int m_Radius;
  unsigned int m_Direction;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDirectionalRectangularSmoothFilter.txx"
#endif

#endif
