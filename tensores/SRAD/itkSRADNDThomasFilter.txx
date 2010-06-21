/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSRADNDThomasFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkSRADNDThomasFilter_txx
#define _itkSRADNDThomasFilter_txx

#include "itkSRADNDThomasFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
SRADNDThomasFilter<TInputImage, TOutputImage>::SRADNDThomasFilter()
{
	m_TimeStep    = 2.0f;
	this->SetNumberOfRequiredInputs( 2 );
}

template < class TInputImage, class TOutputImage >
void SRADNDThomasFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	this->Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	InputImagePointer     inputPtr    =   const_cast< TInputImage * >( this->GetInput()  );
	InputImagePointer     tensorPtr   =   const_cast< TInputImage * >( this->GetTensor() );
	OutputImagePointer    outputPtr   =   this->GetOutput();
	
	if ( !inputPtr || !tensorPtr || !outputPtr )
		return;
	
	// Get a copy of the input requested region (should equal the output requested region)
	InputImageRegionType inputRequestedRegion = inputPtr->GetRequestedRegion();
	
	// Pad the input requested region by the operator radius
	InputSizeType radius;
	radius.Fill( 1 );
	inputRequestedRegion.PadByRadius( radius                                );
	inputRequestedRegion.Crop(        inputPtr->GetLargestPossibleRegion()  );
	inputPtr->SetRequestedRegion(     inputRequestedRegion                  );
	tensorPtr->SetRequestedRegion(    inputRequestedRegion                  );
    return;
}


template < class TInputImage, class TOutputImage >
void SRADNDThomasFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	//-------------------------------------------------------------------------------------
	// Boundary conditions for this filter:
	ZeroFluxNeumannBoundaryCondition<InputImageType>  nbc;	
	// Iterators:
	ConstNeighborhoodIterator<InputImageType>  iit; // Input image
	ConstNeighborhoodIterator<InputImageType>  tit; // Tensor image
	ImageRegionIterator<OutputImageType>       oit; // Output image
	// Allocate input and output
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input   =  this->GetInput();
	InputImageConstPointer   tensor  =  this->GetTensor();
	// Offset object for the input iterator:
	typename ConstNeighborhoodIterator< InputImageType >::OffsetType offset1;
	typename ConstNeighborhoodIterator< InputImageType >::OffsetType offset2;
	
	// Find the data-set boundary "faces"
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType           faceList;
	NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>                                  bC;
	InputSizeType radius;
	radius.Fill( 1 );
	faceList = bC( input, outputRegionForThread, radius );
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
	
	unsigned int DIV = 2*(TInputImage::ImageDimension);
	for (fit=faceList.begin(); fit != faceList.end(); ++fit){ // Iterate through facets
		// Iterators:
		iit   = ConstNeighborhoodIterator<InputImageType>(  radius, input,  *fit  );
		tit   = ConstNeighborhoodIterator<InputImageType>(  radius, tensor, *fit  );
		oit   = ImageRegionIterator<OutputImageType>(       output,         *fit  );
		// Boundary condition:
		iit.OverrideBoundaryCondition( &nbc );
		tit.OverrideBoundaryCondition( &nbc );
		// Iterations:
		for( iit.GoToBegin(), tit.GoToBegin(), oit.GoToBegin(); !iit.IsAtEnd(); ++iit,++tit,++oit ){
			double pixel = 0.0f;
			for( unsigned int k=0; k<TInputImage::ImageDimension; k++ ){
				offset1.Fill( 0 );
				offset2.Fill( 0 );
				offset1[k] = 1;
				offset2[k] = 1;
				pixel += ((double)(tit.GetPixel(offset1)))*(((double)(iit.GetPixel(offset2)))-((double)(iit.GetCenterPixel())));
				offset1[k] = 0;
				offset2[k] = -1;
				pixel += ((double)(tit.GetPixel(offset1)))*(((double)(iit.GetPixel(offset2)))-((double)(iit.GetCenterPixel())));
			}
			pixel = m_TimeStep/((double)DIV)*pixel + ( (double)(iit.GetCenterPixel()) );
			// Set the output pixel:
			oit.Set(   static_cast<OutputPixelType>( pixel )   );			
		}
	}
}

template < class TInputImage, class TOutputImage >
void SRADNDThomasFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "TimeStep:    " << m_TimeStep    << std::endl;
}

} // end namespace itk

#endif
