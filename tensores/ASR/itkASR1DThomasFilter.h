/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkASR1DThomasFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkASR1DThomasFilter_h
#define __itkASR1DThomasFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ASR1DThomasFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef ASR1DThomasFilter                               Self;
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
	itkTypeMacro( ASR1DThomasFilter, ImageToImageFilter );

	/** Types for the diffusion tensor field */
	static const unsigned int TensorDimension = ( (TInputImage::ImageDimension)*(TInputImage::ImageDimension + 1) )/2;
	typedef itk::FixedArray< float, TensorDimension >       TensorType;
	typedef itk::Image< TensorType, TInputImage::ImageDimension > 
		                                                    TensorImageType;
	typedef typename TensorImageType::Pointer               TensorImagePointer;
	typedef typename TensorImageType::ConstPointer          TensorImageConstPointer;

	/** Set and get the diffusion tensor field */
	void SetTensorField( TensorImagePointer tensor ){
		m_TensorField = tensor;
		return;
	}
	TensorImagePointer GetTensorField( void ){
		return( m_TensorField );
	}
	/** Set and get the time constant */
	itkSetMacro(                 TimeStep,   float        );
	itkGetConstReferenceMacro(   TimeStep,   float        );
	/** Set and get the constant value of iso-contour eigenvalue */
	itkSetMacro(                 Alpha,      float        );
	itkGetConstReferenceMacro(   Alpha,      float        );
	/** Set and get the filtering dimension */
	itkSetMacro(                 Direction,  unsigned int );
	itkGetConstReferenceMacro(   Direction,  unsigned int );
protected:
	ASR1DThomasFilter();
	virtual ~ASR1DThomasFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	int  SplitRequestedRegion( int i, int num, OutputImageRegionType& splitRegion );
	void EnlargeOutputRequestedRegion( DataObject *output );
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
	ASR1DThomasFilter(const Self&);         // Purposely not implemented
	void operator=(const Self&);            // Purposely not implemented
	float              m_TimeStep;          // The time constant for anisotropic diffusion
	float              m_Alpha;             // The constant value of iso-contour direction eigenvalue
	TensorImagePointer m_TensorField;       // The diffusion tensor field
	unsigned int       m_Direction;         // The dimension to filter
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkASR1DThomasFilter.txx"
#endif

#endif
