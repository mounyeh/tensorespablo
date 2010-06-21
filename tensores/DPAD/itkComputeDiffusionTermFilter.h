/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeDiffusionTermFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComputeDiffusionTermFilter_h
#define __itkComputeDiffusionTermFilter_h

#include "itkImageToImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkDirectionalRectangularSmoothFilter.h"
#include "itkComputeCoefficientOfVariation.h"
#include "itkCoefficientOfVariationToDiffusionTerm.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ComputeDiffusionTermFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef ComputeDiffusionTermFilter                      Self;
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
	itkTypeMacro( ComputeDiffusionTermFilter, ImageToImageFilter );


	/** Types for the auxiliar filters */
	typedef itk::CastImageFilter< InputImageType, OutputImageType >
															CastType;
	typedef typename CastType::Pointer                      CastPointer;
	typedef itk::SquareImageFilter< OutputImageType, OutputImageType >
		                                                    SquareType;
	typedef typename SquareType::Pointer                    SquarePointer;
	typedef itk::DirectionalRectangularSmoothFilter< OutputImageType, OutputImageType >
		                                                    SmoothType;
	typedef typename SmoothType::Pointer                    SmoothPointer;
	typedef itk::ComputeCoefficientOfVariation< OutputImageType, OutputImageType >
															CoefficientOfVariationType;
	typedef typename CoefficientOfVariationType::Pointer    CoefficientOfVariationPointer;
	typedef itk::CoefficientOfVariationToDiffusionTerm< OutputImageType, OutputImageType >
															DiffusionType;
	typedef typename DiffusionType::Pointer                 DiffusionPointer;

	void SetUseModeDiffusionTerm( DiffusionTermMode mode )
	{
		m_DiffusionTermMode = mode;
	}
	DiffusionTermMode GetUseModeDiffusionTerm( void ){
		return m_DiffusionTermMode;
	}	
	void SetUseAjaDiffusionTerm( void )
	{
		m_DiffusionTermMode = __Aja__;
	}
	void SetUseYuDiffusionTerm( void )
	{
		m_DiffusionTermMode = __Yu__;
	}
	void SetUseSimplifiedDiffusionTerm( void )
	{
		m_DiffusionTermMode = __Simplified__;
	}
	
	/** Set and get the parameters to compute the diffusion Term: */
	itkSetMacro(                 Radius,     InputSizeType  );
	itkGetConstReferenceMacro(   Radius,     InputSizeType  );
	
	itkSetMacro(                 UseMode,    bool  );
	itkGetConstReferenceMacro(   UseMode,    bool  );
	void SetUseModeOn(){
		this->SetUseMode( true );
	}
	void SetUseModeOff(){
		this->SetUseMode( false );
	}
protected:
	ComputeDiffusionTermFilter();
	virtual ~ComputeDiffusionTermFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void EnlargeOutputRequestedRegion(DataObject *output);
	void GenerateData( );
	double FindMode( double* buffer, unsigned long count );
	unsigned long ComputeHistogram( double* data, unsigned long count, double lower_limit, double upper_limit, float* histogram, unsigned int bins );
	void SmoothHistogram( float* histogram_before, float* histogram_after, unsigned int bins, float* filter, unsigned int filter_length );
	void GaussWin( float* filter, unsigned int radius );
	double QuickSort( double* list, unsigned long init, unsigned long end, unsigned long initial_divisor, unsigned long target );
private:
	ComputeDiffusionTermFilter(const Self&);       // Purposely not implemented
	InputSizeType     m_Radius;                    // Vicinity radius (for future use in subclasses)
	bool              m_UseMode;
	DiffusionTermMode m_DiffusionTermMode;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeDiffusionTermFilter.txx"
#endif

#endif
