/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSRADDiffusionTensorFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkSRADDiffusionTensorFilter_txx
#define _itkSRADDiffusionTensorFilter_txx

#include "itkSRADDiffusionTensorFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
SRADDiffusionTensorFilter<TInputImage, TOutputImage>::SRADDiffusionTensorFilter()
{
	m_Beta       = 0.02f;
}

template < class TInputImage, class TOutputImage >
void SRADDiffusionTensorFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	this->Superclass::Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	InputImagePointer     inputPtr    =   const_cast< TInputImage * >( this->GetInput() );
	OutputImagePointer    outputPtr   =   this->GetOutput();
	
	if ( !inputPtr || !outputPtr )
		return;

	InputImagePointer     firstIt     =   const_cast< TInputImage * >( this->GetFirstIteration() );
	if( !firstIt ){
		this->SetFirstIteration( inputPtr );
		firstIt = inputPtr;
	}
	
	// Get a copy of the input requested region (should equal the output requested region)
	InputImageRegionType inputRequestedRegion = inputPtr->GetRequestedRegion();
	
	// Pad the input requested region by the operator radius
	InputSizeType radius;
	radius.Fill( 1 );
	inputRequestedRegion.PadByRadius( radius );
	inputRequestedRegion.Crop(  inputPtr->GetLargestPossibleRegion()  );
	inputPtr->SetRequestedRegion( inputRequestedRegion );
	firstIt->SetRequestedRegion( inputRequestedRegion );
    return;
}

template < class TInputImage, class TOutputImage >
void SRADDiffusionTensorFilter<TInputImage, TOutputImage >
::GenerateData( )
{
	InputImagePointer firstIt = const_cast< TInputImage * >( this->GetFirstIteration() );
	if( !firstIt )
		this->SetFirstIteration( this->GetInput() );
	this->Superclass::Superclass::GenerateData();
}

template < class TInputImage, class TOutputImage >
void SRADDiffusionTensorFilter<TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	//-------------------------------------------------------------------------------------
	// Boundary conditions for this filter:
	ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;	
	// Iterators:
	ConstNeighborhoodIterator<InputImageType> bit;   // Input image
	ImageRegionConstIterator<InputImageType>  first; // First iteration image
	ImageRegionIterator<OutputImageType>      it;    // Output image
	// Offset object for the input iterator:
	typename ConstNeighborhoodIterator< InputImageType >::OffsetType  offset;
	// Allocate input and output
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input   =  this->GetInput();
	InputImageConstPointer   firstIt =  this->GetFirstIteration();
	
	// Find the data-set boundary "faces"
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType           faceList;
	NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>                                  bC;
	InputSizeType radius;
	radius.Fill( 1 );
	faceList = bC( input, outputRegionForThread, radius );
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
	
	unsigned int NLAP = 2*TInputImage::ImageDimension;

	for (fit=faceList.begin(); fit != faceList.end(); ++fit){ // Iterate through facets
		// Iterators:
		bit   = ConstNeighborhoodIterator<InputImageType>(  radius, input, *fit  );
		first = ImageRegionConstIterator<InputImageType>(   firstIt,    *fit      );
		it    = ImageRegionIterator<OutputImageType>(       output,     *fit      );
		// Boundary condition:
		bit.OverrideBoundaryCondition(&nbc);
		// Iterations:
		for( bit.GoToBegin(),first.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it,++first ){
			// Derivative approximations of grad_R and grad_L and for the laplacian
			double gradR[ TInputImage::ImageDimension ];
			double gradL[ TInputImage::ImageDimension ];
			double lapl = 0.0f;
			for( unsigned int k=0; k<TInputImage::ImageDimension; k++ ){
				gradR[k] = gradL[k] = 0.0f;
				offset.Fill( 0 );
				gradR[k] -= bit.GetPixel( offset );
				gradL[k] += bit.GetPixel( offset );
				offset[k] =  1;
				gradR[k] += bit.GetPixel( offset );
				lapl     += bit.GetPixel( offset );
				offset[k] = -1;
				gradL[k] -= bit.GetPixel( offset );
				lapl     += bit.GetPixel( offset );
			}
			double pixel = static_cast<double>(   bit.GetCenterPixel()   );
			lapl -= ( (double)NLAP )*pixel;
			// Compute q:
			double num = -lapl*lapl/pixel/pixel/16.0f;
			double den = ( 1.0f + lapl/pixel/4.0f )*( 1.0f + lapl/pixel/4.0f );
			for( unsigned int k=0; k<TInputImage::ImageDimension; k++ )
				num += ( gradR[k]*gradR[k] + gradL[k]*gradL[k] )/pixel/pixel/4.0f;
			num /= den;
			if( num<1e-4 )
				num = 0.01f;
			else
				num = ::sqrt( num );
			// Compute c(q):
			num = 1.0f + ( num*num - m_Q0*m_Q0 )/(   m_Q0*m_Q0*( 1.0f + m_Q0*m_Q0 )   );
			num = 1.0f/num;
			// Set the correction:
			double beta = (double)( bit.GetCenterPixel() ) - (double)( first.Get() );
			if( beta<0.0f )
				beta = 1.0f + m_Beta*beta;
			else
				beta = 1.0f - m_Beta*beta;
			if( beta<0.0f ){ beta=0.0f; }
			num *= beta;
			// Set the output pixel:
			it.Set(   static_cast<OutputPixelType>( num )   );			
		}
	}
}



template <class TInputImage, class TOutput>
void SRADDiffusionTensorFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Beta: "       << m_Beta << std::endl;
	os << indent << "Q0:   "       << m_Q0   << std::endl;
}

} // end namespace itk

#endif
