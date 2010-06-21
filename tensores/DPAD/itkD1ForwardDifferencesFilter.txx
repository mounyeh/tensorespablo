/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkD1ForwardDifferencesFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkD1ForwardDifferencesFilter_txx
#define _itkD1ForwardDifferencesFilter_txx

#include "itkD1ForwardDifferencesFilter.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"


namespace itk
{

template <class TInputImage, class TOutputImage>
D1ForwardDifferencesFilter<TInputImage, TOutputImage>::D1ForwardDifferencesFilter()
{
	m_TimeStep      = 2.0f;
	m_DiffusionTerm = static_cast< DiffusionImagePointer >( NULL );
	m_Direction     = 0;
}

template < class TInputImage, class TOutputImage >
void D1ForwardDifferencesFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	// Get pointer to the input
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );	
	if ( !inputPtr ){return;}
	// Generate a padding region:
	InputSizeType radius;
	radius.Fill( 0 );
	radius[m_Direction] = 1;
	// Get a copy of the input requested region:
	InputImageRegionType inputRequestedRegion = inputPtr->GetRequestedRegion();
	// Pad the input requested region:
	inputRequestedRegion.PadByRadius( radius );
	// crop the input requested region at the input's largest possible region
	inputRequestedRegion.Crop(   inputPtr->GetLargestPossibleRegion()   );
	// Set the input requested region:
    inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
}


template< class TInputImage, class TOutputImage >
void D1ForwardDifferencesFilter< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Boundary conditions for this filter:
	itk::ZeroFluxNeumannBoundaryCondition<InputImageType>     nbc1;
	itk::ZeroFluxNeumannBoundaryCondition<DiffusionImageType> nbc2;
	// Iterators:
	itk::ConstNeighborhoodIterator<InputImageType>     bit;  // Input image
	itk::ConstNeighborhoodIterator<DiffusionImageType> dit;  // Diffusion term image
	ImageRegionIterator<OutputImageType>               it;   // Output image
	// Allocate input and output
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input   =  this->GetInput();
	
	// Find the data-set boundary "faces"
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType           faceList;
	NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>                                  bC;
	InputSizeType radius;
	radius.Fill( 0 );
	radius[m_Direction] = 1;
	faceList = bC( input, outputRegionForThread, radius );
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
	
	float DIV = 1.0f/(float)TInputImage::ImageDimension;
	
	for (fit=faceList.begin(); fit != faceList.end(); ++fit){ // Iterate through facets
		// Iterators:
		bit = ConstNeighborhoodIterator<InputImageType>(      radius, input,           *fit  );
		dit = ConstNeighborhoodIterator<DiffusionImageType>(  radius, m_DiffusionTerm, *fit  );
		it  = ImageRegionIterator<OutputImageType>(           output,                  *fit  );
		// Boundary condition:
		bit.OverrideBoundaryCondition(&nbc1);
		dit.OverrideBoundaryCondition(&nbc2);
		// Iterations:
		bit.GoToBegin();
		for( bit.GoToBegin(), dit.GoToBegin(), it.GoToBegin(); !it.IsAtEnd(); ++bit, ++dit, ++it ){
			double cn    = ( dit.GetPixel( 1 ) + dit.GetPixel( 2 ) );
			double cp    = ( dit.GetPixel( 0 ) + dit.GetPixel( 1 ) );
			double pixel = (   bit.GetPixel( 2 )  -  bit.GetPixel( 1 )   )*cn;
			pixel       -= (   bit.GetPixel( 1 )  -  bit.GetPixel( 0 )   )*cp;
			pixel        = 0.5f*m_TimeStep*pixel + DIV*( bit.GetCenterPixel() );
			it.Set(   static_cast<OutputPixelType>( pixel )   );
		}
	}
}

template <class TInputImage, class TOutput>
void D1ForwardDifferencesFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "TimeStep:      " << m_TimeStep      << std::endl;
	os << indent << "DiffusionTerm: " << m_DiffusionTerm << std::endl;
}

} // end namespace itk

#endif
