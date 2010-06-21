/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDTIEstimateTensorFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/03/27 17:01:10 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDTIEstimateTensorFilter_h
#define __itkDTIEstimateTensorFilter_h

#include "itkImageToImageFilter.h"
#include "itkOtsuStatistics.h"
#include "itkOtsuThreshold.h"
#include"itkThresholdToMaskImageFilter.h"
#include "itkImage.h"
#include "itkDWImages.h"
#include <vector>
#include "itkVector.h"
#include "itkArray.h"

#ifndef DTI_ESTIMATE_DIM
#define DTI_ESTIMATE_DIM 7
#endif

namespace itk
{
/** \class DTIEstimateTensorFilter
 *
 * Least Squares fitting of DWI data; the filter is prepared to perform Weighted 
 * LS by iteratively recomputing the weigths applied to the problem based on 
 * succesive estimates of the tensor components, as proposed by Salvador et al.;
 * if the number of iterations is set to 0, then LS is used, and the estimation
 * is far faster than if iterations are >0, since, the matrix of LS may be
 * precomputed. Note that the user MUST set the number of DWI volumes, since we
 * DO NOT assume a particular image or pixel type, which is, the input image
 * may be a VectorImage as well as an Image obeject with vectorial pixel type. 
 * We DO INCLUDE the baseline T2 image as an unknown in the estimation problem,
 * like in the paper by Salvador et al. However, the estimated baseline is not
 * included as an output at this moment
 *
 * \sa Image
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT DTIEstimateTensorFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef DTIEstimateTensorFilter            Self;

	/** Convenient typedefs for simplifying declarations. */
	typedef TInputImage                           InputImageType;
	typedef typename InputImageType::Pointer      InputImagePointer;
	typedef typename InputImageType::ConstPointer InputImageConstPointer;
	typedef TOutputImage                          OutputImageType;
	typedef typename OutputImageType::Pointer     OutputImagePointer;
	
	typedef itk::DWImages<typename TInputImage::PixelType, 3 > DWIType;
	typedef typename DWIType::Pointer                 DWIPointer;
	typedef typename DWIType::ConstPointer            DWIConstPointer;

	/** Standard class typedefs. */
	typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
	typedef SmartPointer<Self>                                   Pointer;
	typedef SmartPointer<const Self>                             ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro( DTIEstimateTensorFilter, ImageToImageFilter );
  
	/** Image typedef support. */
	typedef typename InputImageType::PixelType           InputPixelType;
	typedef typename OutputImageType::PixelType          OutputPixelType;
	typedef typename OutputImageType::RegionType         OutputImageRegionType;
	
	/** Typedefs for mask computation */
	typedef itk::Image< float, 3 >                                              FeaturesImageType;
	typedef itk::Image< unsigned char, 3 >                                      MaskImageType;
	typedef FeaturesImageType::Pointer                                          FeaturesImagePointer;
	typedef MaskImageType::Pointer                                              MaskImagePointer;
	typedef itk::OtsuStatistics< InputImageType, FeaturesImageType >            StatisticsFilterType;
	typedef typename StatisticsFilterType::Pointer                              StatisticsFilterPointer;
	typedef itk::OtsuThreshold< FeaturesImageType, FeaturesImageType >          ThresholdFilterType;
	typedef typename ThresholdFilterType::Pointer                               ThresholdFilterPointer;
	typedef itk::ThresholdToMaskImageFilter< FeaturesImageType, MaskImageType > MaskFilterType;
	typedef typename MaskFilterType::Pointer                                    MaskFilterPointer;
	
	typedef itk::Image< float, 3 >                                              T2ImageType;
	typedef typename T2ImageType::Pointer                                       T2ImagePointer;
	
	void ComputeThreshold(void);
	void ComputeMask(void);
	void SetMask(MaskImageType* mask){m_Mask=mask;}
	MaskImagePointer GetMask(void){return m_Mask;}
	
	/** Typedefs for matrix computations; we have a variable number of gradients,
	  * so we do need growing arrays; but gradients are always represented by
	  * means of 3-component vectors; in general, we use fixed-size arrays where
	  * possible, and dynamic size arrays only if needed
	 */
	// Type of the gradients direction:
	typedef itk::Vector<double,3>                    GradientType;
	// Type for each equation line:
	typedef itk::Vector<double,DTI_ESTIMATE_DIM>     EquationType;
	// Type for the storage of the matrix X of the least squares problem:
	typedef std::vector< EquationType >              LSMatrixType;
	// Type for the precomputed inverse of the LS matrix:
	typedef itk::Matrix< double, DTI_ESTIMATE_DIM, DTI_ESTIMATE_DIM >
							 InverseLSType;
	// Type for the vector of unknowns, which is, the components of the tensor:
	typedef itk::Vector< double, DTI_ESTIMATE_DIM >  UnknownsType;
	// baselines may be present in the data):
	typedef itk::Array<unsigned int>                 IndicatorType;
  
	/** Set and get the parameters */
	itkSetMacro( Channels, unsigned int );
	itkGetMacro( Channels, unsigned int );
	itkSetMacro( NBaselines, unsigned int );
	itkGetMacro( NBaselines, unsigned int );
	itkSetMacro( B, double );
	itkGetMacro( B, double );
	itkSetMacro( Iterations, unsigned int );
	itkGetMacro( Iterations, unsigned int );
	itkSetMacro( Threshold, double );
	itkGetMacro( Threshold, double );
	itkSetMacro( ComputeT2, bool );
	itkGetMacro( ComputeT2, bool );
	itkSetMacro( T2, T2ImagePointer );
	itkGetMacro( T2, T2ImagePointer );
	
	/** Add a new gradient direction: */
	void AddGradientDirection( GradientType grad ){
		if( m_B < 0 )
			itkExceptionMacro( << "Please, set the b parameter before!" );
		// Create a vector with the new equation:
		EquationType eq;
		eq[0] =  0.05f;
		eq[1] = -0.1f * grad[0] * grad[0];
		eq[2] = -0.1f * grad[0] * grad[1] * 2;
		eq[3] = -0.1f * grad[0] * grad[2] * 2;
		eq[4] = -0.1f * grad[1] * grad[1];
		eq[5] = -0.1f * grad[1] * grad[2] * 2;
		eq[6] = -0.1f * grad[2] * grad[2];
		// Add the equation to the set:
		m_X.push_back( eq );
		// We do not touch the value of m_Channels; it is the responsibility of
		// the user to set this value properly.
		return;
	}
	
	/** Set the vector with the DWI channels that are going to be used: */
	void SetDWIChannels( IndicatorType ind ){
		m_Indicator = ind;
	}
	void SetBaselines( IndicatorType ind ){
		m_Baselines = ind;
	}
	
	/** Set the vector of DWI channels using c-style vector. The user must set
	 *  m_Channels before */
	void SetDWIChannels( unsigned int* ind ){
		m_Indicator.SetSize( m_Channels );
		for( unsigned int k=0; k<m_Channels; ++k )
			m_Indicator[k] = ind[k];
	}
	void SetBaselines( unsigned int* ind ){
		m_Baselines.SetSize( m_NBaselines );
		for( unsigned int k=0; k<m_NBaselines; ++k )
			m_Baselines[k] = ind[k];
	}
	IndicatorType GetDWIChannels(void){ return m_Indicator; }
	IndicatorType GetBaselines(void){ return m_Baselines; }
	// Custom method to use with itk::DWImages
	void Configure( void );
protected:
	DTIEstimateTensorFilter();
	virtual ~DTIEstimateTensorFilter() {}
	// Threaded filter!
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
	void BeforeThreadedGenerateData();
	// Method to compute the inverse of a symmetric, regular, matrix.
	bool ComputeInverseMatrix( InverseLSType input, InverseLSType& output ) const;
private:
	DTIEstimateTensorFilter(const Self&);  // purposely not implemented
	void operator=(const Self&);           // purposely not implemented
	// The number of DWI channels to use:
	unsigned int m_Channels;
	unsigned int m_NBaselines;
	// The b parameter:
	double       m_B;
	// The number of iterations to refine the WLS weights (0 to use simple LS)
	unsigned int m_Iterations;
	// The matrix where we store X, the matrix of the LS problem:
	LSMatrixType         m_X;
	// The matrix to store the precomputed inverse for the simple LS problem:
	InverseLSType        m_Inverse;
	// The precomputed pseudo-inverse for the LS problem; it should be 7xN, but
	// we store its transpose for convenience:
	LSMatrixType         m_PseudoInverse;
	// The indicator of which channels correspond to useful DWI channels:
	IndicatorType        m_Indicator;
	IndicatorType        m_Baselines;
	// The threshold to compute the mask:
	double               m_Threshold;
	// Mask image:
	MaskImagePointer     m_Mask;
	FeaturesImagePointer m_Features;
	// T2 image:
	T2ImagePointer       m_T2;
	bool                 m_ComputeT2;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDTIEstimateTensorFilter.txx"
#endif

#endif
