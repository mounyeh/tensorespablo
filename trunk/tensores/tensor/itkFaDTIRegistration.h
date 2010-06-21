#ifndef __FaDTIRegistration_h
#define __FaDTIRegistration_h

#include "itkNumericTraitsDTIPixel.h"
#include "itkProcessObject.h"
#include "itkComputeDTIScalars.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkMaskCCBlockMatchingFilter.h"
#include "itkMaskBayesianRegularizationFilter.h"
#include "../registrado_tristan/MLXB3D/itkMAPDisplacementFilter.h"
#include "../registrado_tristan/MLXB3D/itkGaussianSmoothFilter.h"
#include "itkD3BSplineMaskTransform.h"
#include "itkMaskResampleImageFilter.h"
#include "itkResampleImageFilter.h"

//#include "../registrado_tristan/MLXB3D/itkCCBlockMatchingFilter.h"
//#include "../registrado_tristan/MLXB3D/itkBayesianRegularizationFilter.h"
//#include "../registrado_tristan/MLXB3D/itkD3BSplineTransform.h"

#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkModifiedLinearInterpolateImageFunction.h"

namespace itk
{

template < class TInputImage1, class TInputImage2, class TLabeledImage >
class ITK_EXPORT FaDTIRegistration : public ProcessObject
{
	public:
		/** Standard class typedefs. */
		typedef FaDTIRegistration< TInputImage1, TInputImage2, TLabeledImage >			Self;
		typedef ProcessObject                                                           Superclass;
		typedef SmartPointer<Self>                                                      Pointer;	
		typedef SmartPointer<const Self>                                                ConstPointer;
	
		/** Convenient typedefs for simplifying declarations. */
		typedef TInputImage1                                                            FixedImageType;
		typedef typename FixedImageType::Pointer                                        FixedImagePointer;
		typedef TInputImage2                                                            MovingImageType;
		typedef typename MovingImageType::Pointer                                       MovingImagePointer;
		typedef TLabeledImage															LabelImageType;
		typedef typename LabelImageType::Pointer										LabelImagePointer;
		typedef TLabeledImage															OutputImageType;
		typedef typename OutputImageType::Pointer                                       OutputImagePointer;
	
		typedef typename TInputImage1::SizeType                                          SizeType ;
		typedef typename TInputImage2::PointType                                         PointType;

		typedef typename FixedImageType::PixelType										PixelType;
		//typedef typename itk::NumericTraits<PixelType>::ScalarRealType					InternalPixelType;
		typedef float																	InternalPixelType;
		typedef itk::Image<InternalPixelType,3>											InternalImageType;
		typedef typename InternalImageType::Pointer										InternalImagePointer;
		
		typedef itk::Image<bool,3> MaskImageType;
		
		// Compute FA filters
		typedef itk::ComputeDTIScalars<FixedImageType,InternalImageType>				ComputeScalarsType;
		typedef typename ComputeScalarsType::Pointer						    		ComputeScalarsPointer;
 
		// Multiresolution pyramid
		typedef itk::RecursiveMultiResolutionPyramidImageFilter< InternalImageType, InternalImageType >
																						 PyramidType;
		typedef typename PyramidType::Pointer                                            PyramidPointer;

		typedef itk::MultiResolutionPyramidImageFilter< MaskImageType, MaskImageType >
																						 MaskPyramidType;
		typedef typename MaskPyramidType::Pointer                                        MaskPyramidPointer;


		typedef itk::MaskCCBlockMatchingFilter< InternalImageType, InternalImageType, MaskImageType >       CCFilterType;
	//	typedef itk::CCBlockMatchingFilter< InternalImageType, InternalImageType >       CCFilterType;
		typedef typename CCFilterType::Pointer                                           CCFilterPointer; 
		     
		typedef itk::MaskBayesianRegularizationFilter< MaskImageType, FixedImageType::ImageDimension >      RegFilterType;
	//	typedef itk::BayesianRegularizationFilter< FixedImageType::ImageDimension >      RegFilterType;
		typedef typename RegFilterType::Pointer                                          RegFilterPointer;
		
		typedef itk::MAPDisplacementFilter< PointType , FixedImageType::ImageDimension > MAPFilterType;
		typedef typename MAPFilterType::Pointer                                          MAPFilterPointer;	
		
		typedef itk::GaussianSmoothFilter< double, FixedImageType::ImageDimension >      SmoothFilterType;
		typedef typename SmoothFilterType::Pointer                                       SmoothFilterPointer;
		
		typedef itk::D3BSplineMaskTransform< double, FixedImageType::ImageDimension>	TransformType;
	//	typedef itk::D3BSplineTransform< double, FixedImageType::ImageDimension>	TransformType;
		typedef typename TransformType::Pointer                                          TransformPointer;
		typedef typename TransformType::ParametersType                                         ParametersType;

		typedef itk::MaskResampleImageFilter< InternalImageType, InternalImageType, MaskImageType >         ResampleFilterType;
	//	typedef itk::ResampleImageFilter< InternalImageType, InternalImageType >         ResampleFilterType;
		typedef typename ResampleFilterType::Pointer                                     ResampleFilterPointer;

		typedef itk::ModifiedLinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
		typedef typename InterpolatorType::Pointer                                       InterpolatorPointer;
	
		typedef itk::MaskResampleImageFilter< LabelImageType, LabelImageType, MaskImageType >				 FinalResampleType;
	//	typedef itk::ResampleImageFilter< LabelImageType, LabelImageType >				 FinalResampleType;
		typedef typename FinalResampleType::Pointer										 FinalResamplePointer;
	
		typedef itk::NearestNeighborInterpolateImageFunction< LabelImageType, double >			 FinalInterpolatorType;
		typedef typename FinalInterpolatorType::Pointer                                  FinalInterpolatorPointer;
	
	
		/** Method for creation through the object factory. */
		itkNewMacro( Self );

		/** Run-time type information (and related methods). */
		itkTypeMacro( FaDTIRegistration, ProcessObject );
		
		/** Set the first input. */
		void SetFixed( FixedImageType * image ){
			m_Fixed = image;
		}
	
		/** Set the second input. */
		void SetMoving( MovingImageType * image ){
			m_Moving = image;
		}

		/** Set the third input. */
		void SetLabels( LabelImageType * image ){
			m_Labels = image;
		}

		/** Get the first input. */
		FixedImagePointer  GetFixed( void ){
			return m_Fixed;
		}
  
		/** Get the second input. */
		MovingImagePointer GetMoving( void ){
			return m_Moving;
		}

		/** Get the third input. */
		LabelImagePointer GetLabels( void ){
			return m_Labels;
		}

		/** Get the output. */
		OutputImagePointer GetOutput( void ){
			return m_Output;
		}

		/** Set and get the number of registration levels */
		itkSetMacro(               NLevels, unsigned int      );
		itkGetConstReferenceMacro( NLevels, unsigned int      );

		/** Set and get the downsampling rate */
		itkSetMacro(               Sample,     SizeType       );
		itkGetConstReferenceMacro( Sample,     SizeType       );

		/** Set and get the block size */
		itkSetMacro(               BlockSize, SizeType        );
		itkGetConstReferenceMacro( BlockSize, SizeType        );

		/** Set and get the search size */
		itkSetMacro(               SearchSize, SizeType       );
		itkGetConstReferenceMacro( SearchSize, SizeType       );

		/** Set and get the squared variance */	
		itkSetMacro(               Sigma, itk::Array<double>  );
		itkGetConstReferenceMacro( Sigma, itk::Array<double>  );
	
		/** Set and get the number of regularization steps */
		itkSetMacro(               NIter, unsigned int        );
		itkGetConstReferenceMacro( NIter, unsigned int        );

		/** Set and get the width of the gaussian kernel */
		itkSetMacro(               GaussianRadius, SizeType   );
		itkGetConstReferenceMacro( GaussianRadius, SizeType   );
	
		/** Set and get the mixing parameter */
		itkSetMacro(               Lambda,   double           );
		itkGetConstReferenceMacro( Lambda,   double           );

		/** Set and get the final trasformation parameters */	
		itkSetMacro(               Parameters, ParametersType);
		itkGetConstReferenceMacro( Parameters, ParametersType  );
	
		/** Start registration */
		void Start( void );

	protected:
		FaDTIRegistration();
		virtual ~FaDTIRegistration() {}
		void PrintSelf(std::ostream& os, Indent indent) const;
	
	private:
		FaDTIRegistration(const Self&);                // Purposely not implemented
		void operator=(const Self&);                  // Purposely not implemented

		unsigned int            m_NLevels;            // Number of registration levels
		SizeType                m_Sample;             // Down-sampling rate across each direction
		SizeType                m_BlockSize;          // Radius of the block to match
		SizeType                m_SearchSize;         // Radius of the block matching search vicinity
		itk::Array<double>      m_Sigma;              // Variance of the algorithm 
		unsigned int            m_NIter;              // Number of regularization steps
		SizeType                m_GaussianRadius;     // Width of the smoothing gaussian kernel over each dimension
		double                  m_Lambda;             // Mixing proportion between forward and inverse transforms

		FixedImagePointer       m_Fixed;
		MovingImagePointer      m_Moving;
		LabelImagePointer		m_Labels;
		
		OutputImagePointer      m_Output;
		
		typename TransformType::ParametersType m_Parameters;

	
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFaDTIRegistration.txx"
#endif

#endif
