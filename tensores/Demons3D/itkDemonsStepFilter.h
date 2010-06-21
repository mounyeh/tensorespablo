/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsStepFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDemonsStepFilter_h
#define __itkDemonsStepFilter_h

#include "itkProcessObject.h"
#include "itkImage.h"

#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkLocalMeanFilter.h"
#include "itkSmoothFieldFilter.h"
#include "itkComputeLCCFilter.h"
#include "itkComputeCorrectionFilter.h"
#include "itkAddVectorsFilter.h"
#include "itkWarpImageFilter.h"
#include "itkModifiedLinearInterpolateImageFunction.h"
#include "itkSquareImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"

namespace itk
{

template < const unsigned int NDimension >
class ITK_EXPORT DemonsStepFilter : public ImageToImageFilter<  itk::Image< double, NDimension >, itk::Image< itk::Vector< double, NDimension >, NDimension >  >
{
public:
	/** Convenient typedefs for simplifying declarations. */
	typedef itk::Image< double, NDimension >                                         InputImageType;
	typedef typename InputImageType::Pointer                                         InputImagePointer;
	typedef itk::Image< itk::Vector< double, NDimension >, NDimension >              OutputImageType;
	typedef typename OutputImageType::Pointer                                        OutputImagePointer;

	/** Standard class typedefs. */
	typedef DemonsStepFilter< NDimension >                                           Self;
	typedef ImageToImageFilter< InputImageType, OutputImageType >                    Superclass;
	typedef SmartPointer< Self >                                                     Pointer;
	typedef SmartPointer< const Self >                                               ConstPointer;
	
	typedef typename InputImageType::PointType                                       PointType;

	typedef itk::Image< itk::Vector< double, NDimension >, NDimension >              DeformationFieldType;
	typedef typename DeformationFieldType::Pointer                                   DeformationFieldPointer;

	typedef itk::LocalMeanFilter< double, NDimension >                               LocalStatsFilterType;
	typedef itk::SmoothFieldFilter< double, NDimension >                             SmoothingFilterType;
	typedef itk::GradientRecursiveGaussianImageFilter< InputImageType, 
		                                               DeformationFieldType >        GradientCalculatorType;
	typedef itk::ComputeLCCFilter< NDimension >                                      ComputeLCCType;
	typedef itk::ComputeCorrectionFilter< NDimension >                               ComputeCorrectionType;
	typedef itk::WarpImageFilter<   InputImageType, 
		                            InputImageType, 
		                            DeformationFieldType   >                         WarpFilterType;
	typedef itk::ModifiedLinearInterpolateImageFunction< InputImageType, double >    InterpolateType;


	typedef typename LocalStatsFilterType::Pointer                                   LocalStatsFilterPointer;
	typedef typename SmoothingFilterType::Pointer                                    SmoothingFilterPointer;
	typedef typename GradientCalculatorType::Pointer                                 GradientCalculatorPointer;
	typedef typename ComputeLCCType::Pointer                                         ComputeLCCPointer;
	typedef typename ComputeCorrectionType::Pointer                                  ComputeCorrectionPointer;
	typedef typename WarpFilterType::Pointer                                         WarpFilterPointer;
	typedef typename InterpolateType::Pointer                                        InterpolatePointer;

	typedef typename InputImageType::RegionType                                      RegionType ;
	typedef typename InputImageType::SizeType                                        SizeType ;
	typedef typename InputImageType::IndexType                                       IndexType ;

	
	typedef itk::AddVectorsFilter< double, NDimension >                              AddFilterType;
	typedef typename AddFilterType::Pointer                                          AddFilterPointer;

	typedef MultiplyImageFilter< InputImageType, InputImageType, InputImageType >    CrossProductType;
	typedef typename CrossProductType::Pointer                                       CrossProductPointer;

	typedef SquareImageFilter< InputImageType, InputImageType >                      SquareType;
	typedef typename SquareType::Pointer                                             SquarePointer;	

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( DemonsStepFilter, ImageToImageFilter );

	/** Set and get the Levenberg-Marquardt parameter */
	itkSetMacro(                 Lambda,          double    );
	itkGetConstReferenceMacro(   Lambda,          double    );

	/** Set and get the elastic smoothing variances */
	itkSetMacro(                 SigmaElastic,    SizeType  );
	itkGetConstReferenceMacro(   SigmaElastic,    SizeType  );

	/** Set and get the elastic smoothing variances */
	itkSetMacro(                 SigmaFluid,      SizeType  );
	itkGetConstReferenceMacro(   SigmaFluid,      SizeType  );

	/** Set and get the local stats variances */
	itkSetMacro(                 SigmaStats,      double    );
	itkGetConstReferenceMacro(   SigmaStats,      double    );

	/** Set and get the gradient calculation variances */
	itkSetMacro(                 SigmaGradient,   double    );
	itkGetConstReferenceMacro(   SigmaGradient,   double    );

	//------------------------------------------------------------------------------------
	// Set and get the inputs:
	//     - Input 1: The fixed image I
	//     - Input 2: The moving image J
	//     - Input 3: The smoothed value of I, mean( I )
	//     - Input 4: The smoothed value of I^2, mean( I^2 );
	//     - Input 5: The current deformation field

	/** Set the first input */
	void SetInput1( const InputImageType  * image ){ this->SetInput( 0, image ); }
	/** Set the second input */
	void SetInput2( const InputImageType  * image ){ this->SetInput( 1, image ); }
	/** Set the first input */
	void SetInput3( const InputImageType  * image ){ this->SetInput( 2, image ); }
	/** Set the second input */
	void SetInput4( const InputImageType  * image ){ this->SetInput( 3, image ); }
	/** Get the first input */
	const InputImageType  * GetInput1( void ){ return this->GetInput( 0 ); }
	/** Get the second input */
	const InputImageType  * GetInput2( void ){ return this->GetInput( 1 ); }
	/** Get the first input */
	const InputImageType  * GetInput3( void ){ return this->GetInput( 2 ); }
	/** Get the second input */
	const InputImageType  * GetInput4( void ){ return this->GetInput( 3 ); }
	//------------------------------------------------------------------------------------

	/** Use/not use elastic regularization */
	void SetUseElasticRegularization( bool use ){
		m_UseEReg = use;
	}
	bool GetUseElasticRegularization( void ){
		return m_UseEReg;
	}

	/** Use/not use fluid regularization */
	void SetUseFluidRegularization( bool use ){
		m_UseFReg = use;
	}
	bool GetUseFluidRegularization( void ){
		return m_UseFReg;
	}

	/** Set/get the deformation field */
	DeformationFieldPointer GetDeformationField( void ){
		return m_Field;
	}
	void SetDeformationField( DeformationFieldPointer field ){
		m_Field = field;
	}
protected:
	DemonsStepFilter();
	virtual ~DemonsStepFilter() {}
	void PrintSelf(std::ostream& os, Indent indent) const;
	/** Perform an individual registration step */
	void GenerateData( void );
private:
	DemonsStepFilter(const Self&);                // Purposely not implemented
	void operator=(const Self&);                  // Purposely not implemented

	// Levenberg-Marquardt parameter:
	double                     m_Lambda;

	// Variances for gaussian kernels:
	SizeType                   m_SigmaElastic;
	SizeType                   m_SigmaFluid;
	double                     m_SigmaStats;
	double                     m_SigmaGradient;

	// This allows us to use/not use the regularization step:
	bool                       m_UseEReg;
	bool                       m_UseFReg;

	// The deformation field
	DeformationFieldPointer    m_Field;

	//------------------------------------------------------------------------------------------------
	// The various filters:
	WarpFilterPointer          m_Warp;
	ComputeLCCPointer          m_LCC;
	LocalStatsFilterPointer    m_LStats;
	LocalStatsFilterPointer    m_LJ2Stats;
	SquarePointer              m_Square;
	LocalStatsFilterPointer    m_LIJStats;
	CrossProductPointer        m_Cross;
	GradientCalculatorPointer  m_GradientCalculator;
	ComputeCorrectionPointer   m_Correction;
	AddFilterPointer           m_Add;
	SmoothingFilterPointer     m_FluidReg;
	SmoothingFilterPointer     m_ElasticReg;
	//------------------------------------------------------------------------------------------------
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDemonsStepFilter.txx"
#endif

#endif
