/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMLD.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMLD_h
#define __itkMLD_h

#include "itkProcessObject.h"
#include "itkImage.h"

#include "itkDemonsStepFilter.h"
#include "itkLocalMeanFilter.h"
#include "itkCastImageFilter.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkModifiedLinearInterpolateImageFunction.h"
#include "itkVectorResampleImageFilter.h"
#include "itkSquareImageFilter.h"

namespace itk
{

template < class TInputImage1, class TInputImage2 >
class ITK_EXPORT MLD : public ProcessObject
{
public:
	/** Standard class typedefs. */
	typedef MLD< TInputImage1, TInputImage2 >                                        Self;
	typedef ProcessObject                                                            Superclass;
	typedef SmartPointer<Self>                                                       Pointer;
	typedef SmartPointer<const Self>                                                 ConstPointer;

	/** Convenient typedefs for simplifying declarations. */
	typedef TInputImage1                                                             FixedImageType;
	typedef typename FixedImageType::Pointer                                         FixedImagePointer;
	typedef TInputImage2                                                             MovingImageType;
	typedef typename MovingImageType::Pointer                                        MovingImagePointer;
	typedef TInputImage2                                                             OutputImageType;
	typedef typename OutputImageType::Pointer                                        OutputImagePointer;
	typedef typename TInputImage2::PointType                                         PointType;

	typedef double                                                                   InternalPixelType;
	typedef itk::Image< InternalPixelType, FixedImageType::ImageDimension >          InternalImageType;
	typedef typename InternalImageType::Pointer                                      InternalImagePointer;
	
	typedef itk::CastImageFilter< FixedImageType, InternalImageType >                FixedCastFilterType;
	typedef typename FixedCastFilterType::Pointer                                    FixedCastFilterPointer;
	typedef itk::CastImageFilter< MovingImageType, InternalImageType >               MovingCastFilterType;
	typedef typename MovingCastFilterType::Pointer                                   MovingCastFilterPointer;

	typedef itk::RecursiveMultiResolutionPyramidImageFilter< InternalImageType, InternalImageType >
		                                                                             PyramidType;
	typedef typename PyramidType::Pointer                                            PyramidPointer;

	typedef itk::DemonsStepFilter< FixedImageType::ImageDimension >                  DemonsType;
	typedef typename DemonsType::Pointer                                             DemonsPointer;

	typedef itk::LocalMeanFilter< InternalPixelType, TInputImage1::ImageDimension >  LocalStatsFilterType;
	typedef typename LocalStatsFilterType::Pointer                                   LocalStatsFilterPointer;

	typedef itk::Image<   itk::Vector< double, TInputImage1::ImageDimension >, 
		                  TInputImage1::ImageDimension    >                          DeformationFieldType;
	typedef typename DeformationFieldType::Pointer                                   DeformationFieldTypePointer;
	typedef itk::WarpImageFilter< TInputImage1, TInputImage2, DeformationFieldType > WarpFilterType;
	typedef typename WarpFilterType::Pointer                                         WarpFilterPointer;
	typedef itk::ModifiedLinearInterpolateImageFunction< MovingImageType, double >   InterpolatorType;
	typedef typename InterpolatorType::Pointer                                       InterpolatorPointer;
	typedef itk::VectorResampleImageFilter< DeformationFieldType, DeformationFieldType >   
		                                                                             ResampleFilterType;
	typedef typename ResampleFilterType::Pointer                                     ResampleFilterPointer;

	typedef itk::SquareImageFilter< InternalImageType, InternalImageType >           SquareType;
	typedef typename SquareType::Pointer                                             SquarePointer;

	typedef typename TInputImage1::RegionType                                        RegionType ;
	typedef typename TInputImage1::SizeType                                          SizeType ;
	typedef typename TInputImage1::IndexType                                         IndexType ;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( MLD, ProcessObject );

	/** Set and get the number of registration levels */
	itkSetMacro(                 NLevels,         unsigned int      );
	itkGetConstReferenceMacro(   NLevels,         unsigned int      );

	/** Set and get the number of registration lsteps */
	itkSetMacro(                 Steps,           unsigned int      );
	itkGetConstReferenceMacro(   Steps,           unsigned int      );

	/** Set the initial Levenberg-Marquardt paramter: */
	itkSetMacro(                 Lambda,          double            );
	itkGetConstReferenceMacro(   Lambda,          double            );
	
	/** Set the decreasing rate for lambda: */
	itkSetMacro(                 Tau,             double            );
	itkGetConstReferenceMacro(   Tau,             double            );
	
	/** Set and get the elastic smoothing variances */
	itkSetMacro(                 SigmaElastic,    SizeType          );
	itkGetConstReferenceMacro(   SigmaElastic,    SizeType          );

	void SetSigmaElastic( int s ){
		SizeType size;
		size.Fill( s );
		this->SetSigmaElastic( size );
	}

	/** Set and get the fluid smoothing variances */
	itkSetMacro(                 SigmaFluid,      SizeType          );
	itkGetConstReferenceMacro(   SigmaFluid,      SizeType          );

	void SetSigmaFluid( int s ){
		SizeType size;
		size.Fill( s );
		this->SetSigmaFluid( size );
	}

	/** Set and get the local stats variances */
	itkSetMacro(                 SigmaStats,      double            );
	itkGetConstReferenceMacro(   SigmaStats,      double            );

	/** Set and get the gradient calculation variances */
	itkSetMacro(                 SigmaGradient,   double            );
	itkGetConstReferenceMacro(   SigmaGradient,   double            );

	/** Set the first input. */
	void SetInput1( FixedImageType * image ){
		m_Input1 = image;
	}
	
	/** Set the second input. */
	void SetInput2( MovingImageType * image ){
		m_Input2 = image;
	}

	/** Get the first input. */
	FixedImagePointer  GetInput1( void ){
		return m_Input1;
	}
  
	/** Get the second input. */
	MovingImagePointer GetInput2( void ){
		return m_Input2;
	}

	/** Get the output. */
	OutputImagePointer GetOutput( void ){
		return m_Output;
	}

	/** Start registration */
	void Start( void );

	/** Compute local statistics */
	InternalImagePointer ComputeLocalStatistics( InternalImagePointer image );

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
	
	/** Set and get the deformation field */
	void SetDeformationField(  DeformationFieldTypePointer ptr ){
		m_DeformationField = ptr;
	}
	DeformationFieldTypePointer GetDeformationField(void){
		return m_DeformationField;
	}

protected:
	MLD();
	virtual ~MLD() {}
	void PrintSelf(std::ostream& os, Indent indent) const;
private:
	MLD(const Self&);                             // Purposely not implemented
	void operator=(const Self&);                  // Purposely not implemented

	// Levenberg-Marquardt parameter:
	double                  m_Lambda;
	double                  m_Tau;

	// Variances for gaussian kernels:
	SizeType                m_SigmaElastic;
	SizeType                m_SigmaFluid;
	double                  m_SigmaStats;
	double                  m_SigmaGradient;

	// This allows us to use/not use the regularization step:
	bool                    m_UseEReg;
	bool                    m_UseFReg;

	// Number of registration levels:
	unsigned int            m_NLevels;

	// Number of registration steps:
	unsigned int            m_Steps;
	
	OutputImagePointer      m_Output;
	FixedImagePointer       m_Input1;
	MovingImagePointer      m_Input2;
	
	DeformationFieldTypePointer m_DeformationField;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMLD.txx"
#endif

#endif
