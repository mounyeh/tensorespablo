/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSRADNDThomasFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSRADNDThomasFilter_h
#define __itkSRADNDThomasFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT SRADNDThomasFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef SRADNDThomasFilter                              Self;
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
	itkTypeMacro( SRADNDThomasFilter, ImageToImageFilter );

	/** Types for the diffusion tensor field */
	typedef itk::Image< float,TInputImage::ImageDimension > TensorImageType;
	typedef typename TensorImageType::Pointer               TensorImagePointer;
	typedef typename TensorImageType::ConstPointer          TensorImageConstPointer;

	/** Set and get the diffusion tensor field */
	void SetTensorField( const InputImageType *tensor ){
		this->SetInput( 1, tensor );
		return;
	}
	const InputImageType* GetTensor( void ){
		return( this->GetInput(1) );
	}
	/** Set and get the time constant */
	itkSetMacro(                 TimeStep,   float        );
	itkGetConstReferenceMacro(   TimeStep,   float        );
protected:
	SRADNDThomasFilter();
	virtual ~SRADNDThomasFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
	SRADNDThomasFilter(const Self&);        // Purposely not implemented
	void operator=(const Self&);            // Purposely not implemented
	float              m_TimeStep;          // The time constant for anisotropic diffusion
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSRADNDThomasFilter.txx"
#endif

#endif
