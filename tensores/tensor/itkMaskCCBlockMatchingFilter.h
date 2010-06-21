/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMaskCCBlockMatchingFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:46 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMaskCCBlockMatchingFilter_h
#define __itkMaskCCBlockMatchingFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"
#include "itkArray.h"
#include "itkImage.h"

namespace itk {

template<class TInputImage1, class TInputImage2, class TMaskImage>
class ITK_EXPORT MaskCCBlockMatchingFilter 
	: public ImageToImageFilter<   TInputImage1, itk::Image< itk::Array<double>, TInputImage1::ImageDimension >   >
{
public:
  /** Standard Self typedef */
  typedef MaskCCBlockMatchingFilter                            Self;
  typedef ImageToImageFilter<   TInputImage1, itk::Image< itk::Array<double>, TInputImage1::ImageDimension >   > 
	                                                       Superclass;
  typedef SmartPointer<Self>                               Pointer;
  typedef SmartPointer<const Self>                         ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(MaskCCBlockMatchingFilter, ImageToImageFilter);
  
  /** Image related typedefs. */
  typedef itk::Image< itk::Array<double>, TInputImage1::ImageDimension > OutputImageType;
  typedef typename OutputImageType::Pointer                              OutputImagePointer;
  typedef typename OutputImageType::PixelType                            OutputPixelType;
  typedef typename OutputImageType::RegionType                           OutputImageRegionType;
  typedef typename OutputImageType::SizeType                             OutputImageSizeType;
  typedef typename OutputImageType::IndexType                            OutputImageIndexType;

	typedef          TInputImage1                InputImage1Type;
	typedef          TInputImage2                InputImage2Type;
	typedef          TMaskImage					 MaskImageType;

	typedef typename TInputImage1::Pointer       InputImage1Pointer;
	typedef typename TInputImage2::Pointer       InputImage2Pointer;
	typedef typename TMaskImage::Pointer		 MaskImagePointer;

  typedef typename TInputImage1::ConstPointer  InputImage1ConstPointer;
  typedef typename TInputImage2::ConstPointer  InputImage2ConstPointer;

  typedef typename TInputImage1::RegionType    RegionType ;
  typedef typename TInputImage1::SizeType      SizeType ;
  typedef typename TInputImage1::IndexType     IndexType ;

  typedef typename TInputImage2::RegionType    RegionType2 ;
  typedef typename TInputImage2::SizeType      SizeType2 ;
  typedef typename TInputImage2::IndexType     IndexType2 ;

  typedef typename TInputImage1::PixelType     InputImage1PixelType;
  typedef typename TInputImage2::PixelType     InputImage2PixelType;

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

  /** Set and get the search size */
  itkSetMacro(               SearchSize, SizeType  );
  itkGetConstReferenceMacro( SearchSize, SizeType  );

  /** Set and get the downsampling rate */
  itkSetMacro(               Sample,     SizeType  );
  itkGetConstReferenceMacro( Sample,     SizeType  );
  
  
	/** Set the mask input. */
	void SetMask( MaskImageType * image ){
		m_Mask = image;
	}

	/** Get mask input. */
	MaskImagePointer  GetMask( void ) const{
		return m_Mask;
	}
	
protected:
  MaskCCBlockMatchingFilter();
  ~MaskCCBlockMatchingFilter(){};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Multi-threaded filter */
  void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );

  /** The output is smaller than the input */
  void GenerateOutputInformation( void );

  /** We need a larger input than the output */
  virtual void GenerateInputRequestedRegion( ) throw( InvalidRequestedRegionError );

private:
  MaskCCBlockMatchingFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  // The size of the matching block:
  SizeType m_BlockSize;
  // The size of the search vicinity:
  SizeType m_SearchSize;
  // Downsampling factor in each direction:
  SizeType m_Sample;
  MaskImagePointer m_Mask;
} ; // end of class


} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMaskCCBlockMatchingFilter.txx"
#endif

#endif
