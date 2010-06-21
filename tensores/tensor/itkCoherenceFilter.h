/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCoherenceFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:46 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCoherenceFilter_h
#define __itkCoherenceFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"
#include "itkArray.h"
#include "itkImage.h"

namespace itk {

template<class TInputImage1, class TInputImage2, class TMaskImage>
class ITK_EXPORT CoherenceFilter 
	: public ImageToImageFilter<   TInputImage1, TInputImage1  >
{
public:
	/** Standard Self typedef */
	typedef CoherenceFilter                            Self;
	typedef ImageToImageFilter<   TInputImage1, TInputImage1   > 
																												Superclass;
	typedef SmartPointer<Self>                               Pointer;
	typedef SmartPointer<const Self>                         ConstPointer;
  
	/** Method for creation through the object factory. */
	itkNewMacro(Self);  

	/** Runtime information support. */
	itkTypeMacro(CoherenceFilter, ImageToImageFilter);
  
	/** Image related typedefs. */
	typedef TInputImage1												   OutputImageType;
	typedef typename OutputImageType::Pointer                              OutputImagePointer;
	typedef typename OutputImageType::PixelType                            OutputPixelType;
	typedef typename OutputImageType::RegionType                           OutputImageRegionType;
	typedef          TInputImage1                InputImage1Type;
	typedef          TInputImage2                InputImage2Type;
	typedef          TMaskImage					 MaskImageType;

	typedef typename TInputImage1::Pointer       InputImage1Pointer;
	typedef typename TInputImage2::Pointer       InputImage2Pointer;
	typedef typename TMaskImage::Pointer		 MaskImagePointer;

	typedef typename TInputImage1::ConstPointer  InputImage1ConstPointer;
	typedef typename TInputImage2::ConstPointer  InputImage2ConstPointer;

  typedef typename TInputImage1::SizeType      SizeType ;

  /** Set the first input */
  void SetInput1( const InputImage1Type * image ){
	  this->SetInput( 0, image );
  }
  
  /** Set the second input */
  void SetInput2( const InputImage2Type * image ){
	  this->SetInput( 1, image );
  }

  /** Get the first input */
  const InputImage1Type * GetInput1( void ){
	  return this->GetInput( 0 );
  }
  
  /** Get the second input */
  const InputImage2Type * GetInput2( void ){
	  return this->GetInput( 1 );
  }

 
  /** Image related typedefs */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage1::ImageDimension);

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
  CoherenceFilter();
  ~CoherenceFilter(){};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Multi-threaded filter */
  void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );

  /** The output is smaller than the input */
  //void GenerateOutputInformation( void );

  /** We need a larger input than the output */
  //virtual void GenerateInputRequestedRegion( ) throw( InvalidRequestedRegionError );

private:
  CoherenceFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  // The size of the matching block:
  SizeType m_BlockSize;
  
  MaskImagePointer m_Mask;
  
  unsigned int m_Label;
} ; // end of class


} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCoherenceFilter.txx"
#endif

#endif
