/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkASRNDThomasFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkASRNDThomasFilter_h
#define __itkASRNDThomasFilter_h

#include "itkImageToImageFilter.h"
#include "itkASR1DThomasFilter.h"
#include "itkAddImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ASRNDThomasFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef ASRNDThomasFilter                               Self;
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
	itkTypeMacro( ASRNDThomasFilter, ImageToImageFilter );

	/** Types for the diffusion tensor field */
	static const unsigned int TensorDimension = ( (TInputImage::ImageDimension)*(TInputImage::ImageDimension + 1) )/2;
	typedef itk::FixedArray< float, TensorDimension >       TensorType;
	typedef itk::Image< TensorType, TInputImage::ImageDimension > 
		                                                    TensorImageType;
	typedef typename TensorImageType::Pointer               TensorImagePointer;
	typedef typename TensorImageType::ConstPointer          TensorImageConstPointer;

	/** Types for the auxiliar filters */
	typedef itk::ASR1DThomasFilter< TInputImage, TOutputImage >
		                                                    ThomasType;
	typedef typename ThomasType::Pointer                    ThomasPointer;
	typedef itk::AddImageFilter< TOutputImage, TOutputImage, TOutputImage >
		                                                    AddType;
	typedef typename AddType::Pointer                       AddPointer;

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
protected:
	ASRNDThomasFilter();
	virtual ~ASRNDThomasFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void EnlargeOutputRequestedRegion(DataObject *output);
	void GenerateData( );
private:
	ASRNDThomasFilter(const Self&);         // Purposely not implemented
	void operator=(const Self&);            // Purposely not implemented
	float              m_TimeStep;          // The time constant for anisotropic diffusion
	float              m_Alpha;             // The constant value of iso-contour direction eigenvalue
	TensorImagePointer m_TensorField;       // The diffusion tensor field
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkASRNDThomasFilter.txx"
#endif

#endif
