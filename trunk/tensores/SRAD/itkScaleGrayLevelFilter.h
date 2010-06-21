/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScaleGrayLevelFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScaleGrayLevelFilter_h
#define __itkScaleGrayLevelFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ScaleGrayLevelFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef ScaleGrayLevelFilter                            Self;
	typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
	typedef SmartPointer< Self >                            Pointer;
	typedef SmartPointer< const Self >                      ConstPointer;
	
	/** Convenient typedefs for simplifying declarations. */
	typedef TInputImage                                     InputImageType;
	typedef typename InputImageType::Pointer                InputImagePointer;
	typedef typename InputImageType::ConstPointer           InputImageConstPointer;
	typedef TOutputImage                                    OutputImageType;
	typedef typename OutputImageType::Pointer               OutputImagePointer;
	
	/** Image typedef support. */
	typedef typename InputImageType::PixelType              InputPixelType;
	typedef typename OutputImageType::PixelType             OutputPixelType;
	typedef typename InputImageType::RegionType             InputImageRegionType;
	typedef typename OutputImageType::RegionType            OutputImageRegionType;
	typedef typename InputImageType::SizeType               InputSizeType;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro( ScaleGrayLevelFilter, ImageToImageFilter );

	
	/** Set and get the stopping condition parameter */
	itkSetMacro(                 Factor,       float          );
	itkGetConstReferenceMacro(   Factor,       float          );
protected:
	ScaleGrayLevelFilter();
	virtual ~ScaleGrayLevelFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
	ScaleGrayLevelFilter(const Self&);   // Purposely not implemented
	void operator=(const Self&);         // Purposely not implemented
	float  m_Factor;                     // Eigenvalue correction parameter
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScaleGrayLevelFilter.txx"
#endif

#endif
