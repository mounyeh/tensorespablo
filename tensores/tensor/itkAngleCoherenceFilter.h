/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAngleCoherenceFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:46 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAngleCoherenceFilter_h
#define __itkAngleCoherenceFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"
#include "itkArray.h"
#include "itkImage.h"

namespace itk {

template<class TInputImage, class TMaskImage, class TOutputImage>
class ITK_EXPORT AngleCoherenceFilter 
	: public ImageToImageFilter<   TInputImage, TOutputImage  >
{
public:
	/** Standard Self typedef */
	typedef AngleCoherenceFilter                            Self;
	typedef ImageToImageFilter<   TInputImage, TOutputImage   > 
																												Superclass;
	typedef SmartPointer<Self>                               Pointer;
	typedef SmartPointer<const Self>                         ConstPointer;
  
	/** Method for creation through the object factory. */
	itkNewMacro(Self);  

	/** Runtime information support. */
	itkTypeMacro(AngleCoherenceFilter, ImageToImageFilter);
  
	/** Image related typedefs. */
	typedef TOutputImage												   OutputImageType;
	typedef typename OutputImageType::Pointer                              OutputImagePointer;
	typedef typename OutputImageType::PixelType                            OutputPixelType;
	typedef typename OutputImageType::RegionType                           OutputImageRegionType;
	typedef typename OutputImageType::SizeType                             OutputImageSizeType;
	typedef typename OutputImageType::IndexType                            OutputImageIndexType;

	typedef          TInputImage                InputImageType;
	typedef          TMaskImage					 MaskImageType;

	typedef typename TInputImage::Pointer       InputImagePointer;
	typedef typename TMaskImage::Pointer		 MaskImagePointer;

	typedef typename TInputImage::ConstPointer  InputImageConstPointer;
	typedef typename TInputImage::SizeType      SizeType ;
  
	/** Image related typedefs */
	itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

	/** Set and get the block size */
	itkSetMacro(               BlockSize, SizeType  );
	itkGetConstReferenceMacro( BlockSize, SizeType  );
  
	/** Set the mask input. */
	void SetMask( MaskImageType * image ){
		m_Mask = image;
	}

	/** Get mask input. */
	MaskImagePointer  GetMask( void ) const{
		return m_Mask;
	}

	/** Set and get the label of the region of interest */
	itkSetMacro(               Label, unsigned int      );
	itkGetConstReferenceMacro( Label, unsigned int      );


protected:
  AngleCoherenceFilter();
  ~AngleCoherenceFilter(){};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Multi-threaded filter */
  void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );

private:
  AngleCoherenceFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  // The size of the matching block:
  SizeType m_BlockSize;
  
  MaskImagePointer m_Mask;
  
  unsigned int m_Label;
} ; // end of class


} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAngleCoherenceFilter.txx"
#endif

#endif
