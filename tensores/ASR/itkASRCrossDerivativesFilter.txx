/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkASRCrossDerivativesFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkASRCrossDerivativesFilter_txx
#define _itkASRCrossDerivativesFilter_txx
#include "itkASRCrossDerivativesFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ASRCrossDerivativesFilter<TInputImage, TOutputImage>::ASRCrossDerivativesFilter()
{
	m_TimeStep    = 2.0;
	m_TensorField = static_cast< TensorImagePointer >( NULL );
}

template < class TInputImage, class TOutputImage >
void ASRCrossDerivativesFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();

	// Get pointers to the input and output
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );
	OutputImagePointer outputPtr = this->GetOutput();
	
	if ( !inputPtr || !outputPtr )
		return;
	
	// Get a copy of the input requested region
	InputImageRegionType inputRequestedRegion = inputPtr->GetRequestedRegion();
	// Get a copy of the corresponding region in the tensor field image	
	InputImageRegionType tensorRegion         = m_TensorField->GetBufferedRegion();

	// Pad the region by the 3x3x...x3 neighbourhood
	InputSizeType radius;
	radius.Fill( 1 );
	inputRequestedRegion.PadByRadius(  radius   );

	// Crop the region to their largest possible
	inputRequestedRegion.Crop(      inputPtr->GetLargestPossibleRegion()        );

	// Make sure the tensor field is large enough:
	for( unsigned int k=0; k<TInputImage::ImageDimension; k++ ){
		if( tensorRegion.GetIndex()[k] > inputRequestedRegion.GetIndex()[k] )
			itkExceptionMacro( << "Input requested region outside the buffered region of the tensor field" );
		if( tensorRegion.GetIndex()[k]+tensorRegion.GetSize()[k] < inputRequestedRegion.GetIndex()[k]+inputRequestedRegion.GetSize()[k] )
			itkExceptionMacro( << "Input requested region outside the buffered region of the tensor field" );
	}

	// Set the final requested region
    inputPtr->SetRequestedRegion( inputRequestedRegion );
}


template< class TInputImage, class TOutputImage>
void ASRCrossDerivativesFilter< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	ZeroFluxNeumannBoundaryCondition<InputImageType>  nbc1;
	ZeroFluxNeumannBoundaryCondition<TensorImageType> nbc2;
	
	// Iterators for each image:
	ConstNeighborhoodIterator< InputImageType >     bit;  // Input
	ConstNeighborhoodIterator< TensorImageType >    tit;  // Tensor field
	ImageRegionIterator< OutputImageType >          it;   // Output

	// Offset:
	typename ConstNeighborhoodIterator< InputImageType >::OffsetType  offseti;
	typename ConstNeighborhoodIterator< TensorImageType >::OffsetType offsett;
	
	// Get pointers to the images:
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input   =  this->GetInput();
	
	// Find the data-set boundary "faces"
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator   fit;
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType >::FaceListType           faceList;
	NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType >                                  bC;
	InputSizeType radius;
	radius.Fill( 1 );
	faceList   =   bC( input, outputRegionForThread,radius );
	
	for ( fit=faceList.begin(); fit != faceList.end(); ++fit ){
		// Create iterators for this facet:
		bit = ConstNeighborhoodIterator< InputImageType >(    radius,   input,           *fit   );
		tit = ConstNeighborhoodIterator< TensorImageType >(   radius,   m_TensorField,   *fit   );
		it  = ImageRegionIterator< OutputImageType >(         output,                    *fit   );
		
		// Boundary condition:
		bit.OverrideBoundaryCondition(&nbc1);
		tit.OverrideBoundaryCondition(&nbc2);
		
		// Initiallise:
		bit.GoToBegin();
		tit.GoToBegin();
		
		// For each pixel in this facet:
		while ( ! bit.IsAtEnd() ){
			// Auxiliar output pixel value:
			double value = 0.0;
			// Auxiliar value
			unsigned int pos = 0;
			// Loop to search for the cross-derivatives:
			for( unsigned int f=0; f<TInputImage::ImageDimension-1; f++ ){
				pos++;
				for( unsigned int c=f+1; c<TInputImage::ImageDimension; c++ ){
					//---------------------------------------------------------
					//---------------------------------------------------------
					// Clear offsets
					offseti.Fill( 0 );
					offsett.Fill( 0 );
					// Compute divergence term d_{x_f}( c·d_{x_c}I )
					//---------------------------------------------------------
					offsett[c] =  1; // First term of the divergence derivative
					offseti[c] =  1;
					//---
					offseti[f] =  1; // First term of the gradient derivative
					value += (   (double)( bit.GetPixel(offseti) )   )*(   (double)( tit.GetPixel(offsett)[pos] )   );
					//---
					offseti[f] = -1; // Second term of the gradient derivative
					value -= (   (double)( bit.GetPixel(offseti) )   )*(   (double)( tit.GetPixel(offsett)[pos] )   );
					//---------------------------------------------------------
					offsett[c] = -1; // Second term of the divergence derivative
					offseti[c] = -1;
					//---
					offseti[f] =  1; // First term of the gradient derivative
					value -= (   (double)( bit.GetPixel(offseti) )   )*(   (double)( tit.GetPixel(offsett)[pos] )   );
					//---
					offseti[f] = -1; // Second term of the gradient derivative
					value += (   (double)( bit.GetPixel(offseti) )   )*(   (double)( tit.GetPixel(offsett)[pos] )   );
					//---------------------------------------------------------
					//---------------------------------------------------------
					// Clear offsets
					offseti.Fill( 0 );
					offsett.Fill( 0 );
					// Compute divergence term d_{x_c}( c·d_{x_f}I )
					//---------------------------------------------------------
					offsett[f] =  1; // First term of the divergence derivative
					offseti[f] =  1;
					//---
					offseti[c] =  1; // First term of the gradient derivative
					value += (   (double)( bit.GetPixel(offseti) )   )*(   (double)( tit.GetPixel(offsett)[pos] )   );
					//---
					offseti[c] = -1; // Second term of the gradient derivative
					value -= (   (double)( bit.GetPixel(offseti) )   )*(   (double)( tit.GetPixel(offsett)[pos] )   );
					//---------------------------------------------------------
					offsett[f] = -1; // Second term of the divergence derivative
					offseti[f] = -1;
					//---
					offseti[c] =  1; // First term of the gradient derivative
					value -= (   (double)( bit.GetPixel(offseti) )   )*(   (double)( tit.GetPixel(offsett)[pos] )   );
					//---
					offseti[c] = -1; // Second term of the gradient derivative
					value += (   (double)( bit.GetPixel(offseti) )   )*(   (double)( tit.GetPixel(offsett)[pos] )   );
					//---------------------------------------------------------
					//---------------------------------------------------------
					pos++;
				}
			}
			// Compute the final value:
			value = static_cast< double >(   bit.GetCenterPixel()   ) + (m_TimeStep/4.0f)*value;
			// Set the output pixel value:
			it.Set(    static_cast< OutputPixelType >( value )    );
			// Next pixel:
			++bit; ++tit; ++it;
		}
	}
}

template <class TInputImage, class TOutput>
void ASRCrossDerivativesFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "TimeStep:    " << m_TimeStep    << std::endl;
	os << indent << "TensorField: " << m_TensorField << std::endl;
}

} // end namespace itk

#endif
