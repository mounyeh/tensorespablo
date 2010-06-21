/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianBlurFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkGaussianBlurFilter_txx
#define _itkGaussianBlurFilter_txx
#include "itkGaussianBlurFilter.h"

#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"

#include "math.h"

namespace itk
{

template < class TInputImage, class TOutputImage >
GaussianBlurFilter<  TInputImage, TOutputImage  >
::GaussianBlurFilter()
{
	m_Radius     = 2;
	m_Direction  = 0;
	m_Tol        = 6.0;
	m_Operator.SetSize( (unsigned long)(2*m_Radius+1) );
	m_Correction = 1.6740;
}

template < class TInputImage, class TOutputImage >
void GaussianBlurFilter<  TInputImage, TOutputImage  >
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	InputImagePointer   input  = const_cast< InputImageType* >( this->GetInput() );
	OutputImagePointer  output = this->GetOutput();
	
	if (  (!input) || (!output)  ){ return; }

	// Get a copy of the input requested region (should equal the output requested region)
	InputRegionType inputRequestedRegion = input->GetRequestedRegion();

	SizeType radius;
	radius.Fill( 0 );
	if( m_Direction<TInputImage::ImageDimension )
		radius[m_Direction] = m_Radius;
	else
		itkExceptionMacro( << "Filtering dimension doesn't match with image dimension" );

	inputRequestedRegion.PadByRadius( radius );
	
	if (   inputRequestedRegion.Crop( input->GetLargestPossibleRegion() )   )
		input->SetRequestedRegion( inputRequestedRegion );
	else
		itkExceptionMacro( << "Cannot crop input region to largest possible region" );

	return;
}

template < class TInputImage, class TOutputImage >
void GaussianBlurFilter<  TInputImage, TOutputImage  >
::BeforeThreadedGenerateData( void )
{
	if( m_Radius != 0 ){
		double N   = (double)(2*m_Radius+1);
		double sum = 0.0;
		m_Operator.SetSize( 2*m_Radius+1 );
		m_Correction = 0.0;
		for( unsigned int k=0; k<2*m_Radius+1; k++ ){
			m_Operator[k] = exp(    -(  2.5*( (double)k - (N/2) + 0.5 )/(N/2)  )*(  2.5*( (double)k - (N/2) + 0.5 )/(N/2)  )/2.0    );
		}
		for( unsigned int k=0; k<2*m_Radius+1; k++ ){
			sum += m_Operator[k];
		}
		for( unsigned int k=0; k<2*m_Radius+1; k++ ){
			m_Operator[k] /= sum;
			m_Correction += m_Operator[k];
		}
		m_Correction -= m_Operator[m_Radius];
		m_Correction  = 1/m_Correction;
	}
	else{
		m_Operator.SetSize( 1 );
		m_Operator[0] = 1.0;
		m_Correction  = 1.0;
	}
}

template < class TInputImage, class TOutputImage >
void GaussianBlurFilter<  TInputImage, TOutputImage  >
::ThreadedGenerateData( const OutputRegionType& outputRegionForThread, int threadId )
{
	ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;
	
	// Iterators:
	ConstNeighborhoodIterator<InputImageType>  bit;  // Input
	ImageRegionIterator<OutputImageType>       it;   // Output
	
	// Allocate output
	typename  OutputImageType::Pointer      output = this->GetOutput();
	typename  InputImageType::ConstPointer  input  = this->GetInput();
	
	// Find the data-set boundary "faces"
	SizeType radius;
	radius.Fill( 0 );
	if( m_Direction<TInputImage::ImageDimension )
		radius[m_Direction] = m_Radius;
	else
		itkExceptionMacro( << "Filtering dimension doesn't match with image dimension" );

	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
	NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>                        bC;
	faceList = bC(input, outputRegionForThread, radius);
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
	
	InputPixelType pixel;
	for (fit=faceList.begin(); fit != faceList.end(); ++fit){
		// Iterators:
		bit = ConstNeighborhoodIterator<InputImageType>(  radius, input,  *fit  );
		it  = ImageRegionIterator<OutputImageType>(               output, *fit  );
		//Zero flux boundary condition:
		bit.OverrideBoundaryCondition( &nbc );
		// Initialize pointers:
		it.GoToBegin();
		bit.GoToBegin();
		while ( ! bit.IsAtEnd() ){
			if( m_Radius != 0 ){
				// Causal filtering:
				for( unsigned int d=0; d<TInputImage::ImageDimension; d++ )
					pixel[d] = ( bit.GetPixel(0)[d] )*( m_Operator[0] );
				for( unsigned int k=1; k<(bit.Size()-1)/2; k++ ){
					for( unsigned int d=0; d<TInputImage::ImageDimension; d++ )
						pixel[d] += ( bit.GetPixel(k)[d] )*( m_Operator[k] );
				}
				// Anti-causal filtering:
				for( unsigned int k=(bit.Size()+1)/2; k<bit.Size(); k++ ){
					for( unsigned int d=0; d<TInputImage::ImageDimension; d++ )
						pixel[d] += ( bit.GetPixel(k)[d] )*( m_Operator[k] );
				}
				// Auxiliar values:
				InputPixelType mean;
				for( unsigned int d=0; d<TInputImage::ImageDimension; d++ )
					mean[d] = (pixel[d])*m_Correction;
				InputPixelType center = bit.GetCenterPixel();
				// Outlier rejection:
				for( unsigned int d=0; d<TInputImage::ImageDimension; d++ ){					
					if(   ( mean[d]-center[d] > m_Tol )   ||   ( center[d]-mean[d] > m_Tol )   )
						pixel[d]  = mean[d];
					else
						pixel[d] += ( center[d] )*( m_Operator[m_Radius] );
				}
				// Set the new value:
				it.Set(    static_cast<OutputPixelType>(  pixel            )    );
			}
			else{
				it.Set(    static_cast<OutputPixelType>(  bit.GetPixel(0)  )    );
			}
			
			// Update iterators:
			++bit;
			++it;
		}
	}
}


/**
 * Standard "PrintSelf" method
 */
template < class TInputImage, class TOutputImage >
void GaussianBlurFilter<  TInputImage, TOutputImage  >::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif
