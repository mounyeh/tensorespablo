/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeDiffusionTensorFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComputeDiffusionTensorFilter_h
#define __itkComputeDiffusionTensorFilter_h

#include "itkImageToImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkPartialsToTensorFilter.h"
#include "itkEigenValueCorrectionFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ComputeDiffusionTensorFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef ComputeDiffusionTensorFilter                    Self;
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
	itkTypeMacro( ComputeDiffusionTensorFilter, ImageToImageFilter );

	/** Type for the gradient image */
	typedef itk::CovariantVector< float, TInputImage::ImageDimension > GradientType;
	typedef itk::Image< GradientType, TInputImage::ImageDimension >    GradientImageType;

	/** Types for the auxiliar filters */
	typedef itk::GradientRecursiveGaussianImageFilter< TInputImage, GradientImageType >
		                                                    GradientFilterType;
	typedef typename GradientFilterType::Pointer            GradientFilterPointer;
	typedef itk::PartialsToTensorFilter< GradientImageType, TOutputImage >
		                                                    TensorFilterType;
	typedef typename TensorFilterType::Pointer              TensorFilterPointer;
	typedef itk::EigenValueCorrectionFilter<   TInputImage,   itk::Image< float, TInputImage::ImageDimension >   >
		                                                    CorrectionType;
	typedef typename CorrectionType::Pointer                CorrectionPointer;

	/** Set the first iteration */
	void SetFirstIteration( const InputImageType* firstIteration ){
		this->SetInput( 1, firstIteration );
	}
	const InputImageType* GetFirstIteration( void ){
		return this->GetInput( 1 );
	}

	/** Set and get the parameters to compute the diffusion tensor: */
	itkSetMacro(                 Radius,     InputSizeType  );
	itkGetConstReferenceMacro(   Radius,     InputSizeType  );
	virtual void SetRadius( unsigned int radius ){}
	/** Set and get the width of the gaussian kernels */
	itkSetMacro(                 Sigma,      float          );
	itkGetConstReferenceMacro(   Sigma,      float          );
	/** Set and get the maximum difference between eigenvalues*/
	itkSetMacro(                 Difference, float          );
	itkGetConstReferenceMacro(   Difference, float          );
	/** Set and get the stopping condition parameter */
	itkSetMacro(                 Beta,       float          );
	itkGetConstReferenceMacro(   Beta,       float          );

protected:
	ComputeDiffusionTensorFilter();
	virtual ~ComputeDiffusionTensorFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void EnlargeOutputRequestedRegion(DataObject *output);
	void GenerateData( );
private:
	ComputeDiffusionTensorFilter(const Self&);   // Purposely not implemented
	void operator=(const Self&);                 // Purposely not implemented
	float m_Sigma;                               // The width of the gaussian kernels
	float m_Difference;                          // Parameter to compute the diffusion tensor
	float m_Beta;                                // Eigenvalue correction parameter
	InputSizeType m_Radius;                      // Vicinity radius (for future use in subclasses)
	
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeDiffusionTensorFilter.txx"
#endif

#endif
