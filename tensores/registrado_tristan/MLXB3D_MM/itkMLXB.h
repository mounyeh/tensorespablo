/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMLXB.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMLXB_h
#define __itkMLXB_h

#include "itkProcessObject.h"
#include "itkImage.h"

#include "itkCCBlockMatchingFilter.h"
#include "itkBayesianRegularizationFilter.h"
#include "itkMAPDisplacementFilter.h"
#include "itkGaussianSmoothFilter.h"
#include "itkTransform.h"
#include "itkD3BSplineTransform.h"

#include "itkResampleImageFilter.h"
#include "itkModifiedLinearInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkCastImageFilter.h"

#include "itkRecursiveMultiResolutionPyramidImageFilter.h"

namespace itk
{

template < class TInputImage1, class TInputImage2 >
class ITK_EXPORT MLXB : public ProcessObject
{
public:
	/** Standard class typedefs. */
	typedef MLXB< TInputImage1, TInputImage2 >                                       Self;
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

	typedef float                                                                    InternalPixelType;
	typedef itk::Image< InternalPixelType, FixedImageType::ImageDimension >          InternalImageType;
	
	typedef itk::CastImageFilter< FixedImageType, InternalImageType >                FixedCastFilterType;
	typedef typename FixedCastFilterType::Pointer                                    FixedCastFilterPointer;
	typedef itk::CastImageFilter< MovingImageType, InternalImageType >               MovingCastFilterType;
	typedef typename MovingCastFilterType::Pointer                                   MovingCastFilterPointer;

	typedef itk::RecursiveMultiResolutionPyramidImageFilter< InternalImageType, InternalImageType >
		                                                                             PyramidType;
	typedef typename PyramidType::Pointer                                            PyramidPointer;

	typedef itk::CCBlockMatchingFilter< InternalImageType, InternalImageType >       CCFilterType;
	typedef typename CCFilterType::Pointer                                           CCFilterPointer;      
	typedef itk::BayesianRegularizationFilter< FixedImageType::ImageDimension >      RegFilterType;
	typedef typename RegFilterType::Pointer                                          RegFilterPointer;
	typedef itk::MAPDisplacementFilter< PointType , FixedImageType::ImageDimension > MAPFilterType;
	typedef typename MAPFilterType::Pointer                                          MAPFilterPointer;
	typedef itk::GaussianSmoothFilter< double, FixedImageType::ImageDimension >      SmoothFilterType;
	typedef typename SmoothFilterType::Pointer                                       SmoothFilterPointer;
	typedef itk::D3BSplineTransform< double, FixedImageType::ImageDimension >        TransformType;
	typedef typename TransformType::Pointer                                          TransformPointer;

	typedef itk::ResampleImageFilter< InternalImageType, InternalImageType >         ResampleFilterType;
	typedef typename ResampleFilterType::Pointer                                     ResampleFilterPointer;

	typedef itk::ModifiedLinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
	typedef typename InterpolatorType::Pointer                                       InterpolatorPointer;
	typedef itk::InterpolateImageFunction< InternalImageType, double >               BaseInterpolator;

	typedef itk::ResampleImageFilter< MovingImageType, MovingImageType >             ResampleType;
	typedef typename ResampleType::Pointer                                           ResamplePointer;
	
	typedef itk::LinearInterpolateImageFunction< MovingImageType, double >           FinalInterpolatorType;
	typedef typename FinalInterpolatorType::Pointer                                  FinalInterpolatorPointer;
	typedef itk::InterpolateImageFunction< MovingImageType, double >                 FinalBaseInterpolator;


	

	typedef typename TInputImage1::RegionType                                        RegionType ;
	typedef typename TInputImage1::SizeType                                          SizeType ;
	typedef typename TInputImage1::IndexType                                         IndexType ;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( MLXB, ProcessObject );

	/** Set and get the number of registration levels */
	itkSetMacro(               NLevels, unsigned int      );
	itkGetConstReferenceMacro( NLevels, unsigned int      );

	/** Set and get the block size */
	itkSetMacro(               BlockSize, SizeType        );
	itkGetConstReferenceMacro( BlockSize, SizeType        );
	
	/** Set and get the search size */
	itkSetMacro(               SearchSize, SizeType       );
	itkGetConstReferenceMacro( SearchSize, SizeType       );
	
	/** Set and get the downsampling rate */
	itkSetMacro(               Sample,     SizeType       );
	itkGetConstReferenceMacro( Sample,     SizeType       );

	/** Set and get the squared variance */
	itkSetMacro(               Sigma, itk::Array<double>  );
	itkGetConstReferenceMacro( Sigma, itk::Array<double>  );

	/** Set and get the number of regularization steps */
	itkSetMacro(               NIter, unsigned int        );
	itkGetConstReferenceMacro( NIter, unsigned int        );

	/** Set and get the width of the gaussian kernel */
	itkSetMacro(               GaussianRadius, SizeType   );
	itkGetConstReferenceMacro( GaussianRadius, SizeType   );

	/** Set and get the metric to use */
	itkSetMacro(               Metric, unsigned int       );
	itkGetConstReferenceMacro( Metric, unsigned int       );

	/** Set and get the number of histogram bins for certain metrics */
	itkSetMacro(               Bins, unsigned int         );
	itkGetConstReferenceMacro( Bins, unsigned int         );

	/** Set and get the mixing parameter */
	itkSetMacro(               Lambda,   double           );
	itkGetConstReferenceMacro( Lambda,   double           );

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
protected:
	MLXB();
	virtual ~MLXB() {}
	void PrintSelf(std::ostream& os, Indent indent) const;
private:
	MLXB(const Self&);                            // Purposely not implemented
	void operator=(const Self&);                  // Purposely not implemented

	unsigned int            m_NLevels;            // Number of registration levels
	SizeType                m_BlockSize;          // Radius of the block to match
	SizeType                m_SearchSize;         // Radius of the block matching search vicinity
	SizeType                m_Sample;             // Down-sampling rate across each direction
	itk::Array<double>      m_Sigma;              // Variance of the algorithm 
	SizeType                m_GaussianRadius;     // Width of the smoothing gaussian kernel over each dimension
	unsigned int            m_NIter;              // Number of regularization steps
	unsigned int            m_Metric;             // The metric to use
	unsigned int            m_Bins;               // The number of histogram bins to use
	double                  m_Lambda;             // Mixing proportion between forward and inverse transforms
	
	OutputImagePointer      m_Output;
	FixedImagePointer       m_Input1;
	MovingImagePointer      m_Input2;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMLXB.txx"
#endif

#endif
