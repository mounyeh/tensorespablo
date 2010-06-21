/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkASRFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkASRFilter_h
#define __itkASRFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkASRCrossDerivativesFilter.h"
#include "itkASRNDThomasFilter.h"
#include "itkComputeDiffusionTensorFilter.h"


namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ASRFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef ASRFilter                                       Self;
	typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
	typedef SmartPointer< Self >                            Pointer;
	typedef SmartPointer< const Self >                      ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );
	/** Run-time type information (and related methods). */
	itkTypeMacro( ASRFilter, ImageToImageFilter );
	
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

	/** Internal Image Type */
	typedef itk::Image< float,TInputImage::ImageDimension > InternalImageType;

	/** Types for the diffusion tensor field */
	static const unsigned int TensorDimension = ( (TInputImage::ImageDimension)*(TInputImage::ImageDimension + 1) )/2;
	typedef itk::FixedArray< float, TensorDimension >       TensorType;
	typedef itk::Image< TensorType, TInputImage::ImageDimension > 
		                                                    TensorImageType;
	typedef typename TensorImageType::Pointer               TensorImagePointer;
	typedef typename TensorImageType::ConstPointer          TensorImageConstPointer;
	
	/** Types for the auxiliar filters */
	typedef itk::ASRCrossDerivativesFilter< InputImageType, InternalImageType >    CRType;
	typedef typename CRType::Pointer                                               CRPointer;
	typedef itk::ASRNDThomasFilter< InternalImageType, OutputImageType >           ThomasType;
	typedef typename ThomasType::Pointer                                           ThomasPointer;
	typedef itk::ComputeDiffusionTensorFilter< InputImageType, TensorImageType >   TensorFilterType;
	typedef typename TensorFilterType::Pointer                                     TensorFilterPointer;

	/** Set and get the number of iterations */
	itkSetMacro(               Iter, unsigned int );
	itkGetConstReferenceMacro( Iter, unsigned int );

	/** Set and get the time constant */
	itkSetMacro(               TimeStep, float );
	itkGetConstReferenceMacro( TimeStep, float );
	
	/** Set and get the width of the gaussian kernels */
	itkSetMacro(                 Sigma,      float          );
	itkGetConstReferenceMacro(   Sigma,      float          );
	/** Set and get the maximum difference between eigenvalues*/
	itkSetMacro(                 Difference, float          );
	itkGetConstReferenceMacro(   Difference, float          );
	/** Set and get the stopping condition parameter */
	itkSetMacro(                 Beta,       float          );
	itkGetConstReferenceMacro(   Beta,       float          );
	
	/** Set and get the diffusion tensor filter */
	void SetTensorFilter( TensorFilterType* tensor ){
		m_Tensor = const_cast< TensorFilterType* >( tensor );
	}
	TensorFilterType* GetTensorFilter( void ){
		return m_Tensor;
	}
protected:
	ASRFilter();
	~ASRFilter();
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void EnlargeOutputRequestedRegion( DataObject *output );
	void GenerateData( );
private:
	ASRFilter(const Self&);         // Purposely not implemented
	void operator=(const Self&);    // Purposely not implemented
	unsigned int        m_Iter;     // The number of iterations to perform
	float               m_TimeStep; // Time step for the diffusion process
	TensorFilterPointer m_Tensor;   // The filter to compute the diffusion tensor
	float m_Sigma;                  // The width of the gaussian kernels
	float m_Difference;             // Parameter to compute the diffusion tensor
	float m_Beta;                   // Eigenvalue correction parameter
	InputSizeType m_Radius;         // Vicinity radius (for future use in subclasses)
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkASRFilter.txx"
#endif

#endif

