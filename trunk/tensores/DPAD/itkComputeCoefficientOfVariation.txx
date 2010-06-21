/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeCoefficientOfVariation.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkComputeCoefficientOfVariation_txx
#define _itkComputeCoefficientOfVariation_txx

#include "itkComputeCoefficientOfVariation.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ComputeCoefficientOfVariation<TInputImage, TOutputImage>
::ComputeCoefficientOfVariation()
{
	this->SetNumberOfRequiredInputs( 2 );
	m_NumberOfPixels = 25;
}

template <class TInputImage, class TOutputImage>
void ComputeCoefficientOfVariation<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	// Get pointer to the input
	InputImagePointer  inputPtr1  = const_cast< TInputImage * >( this->GetInput1() );
	InputImagePointer  inputPtr2  = const_cast< TInputImage * >( this->GetInput2() );
	if ( !inputPtr1 || !inputPtr2 ){ return; }
	// This filter needs the entire input:
    inputPtr1->SetRequestedRegion(   inputPtr1->GetRequestedRegion()   );
    inputPtr2->SetRequestedRegion(   inputPtr1->GetRequestedRegion()   );
}

template< class TInputImage, class TOutputImage>
void ComputeCoefficientOfVariation< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Allocate input and output
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input1  =  this->GetInput1();
	InputImageConstPointer   input2  =  this->GetInput2();
	// Iterators:
	itk::ImageRegionConstIterator<InputImageType> ii1 = itk::ImageRegionConstIterator<InputImageType>( input1, outputRegionForThread );
	itk::ImageRegionConstIterator<InputImageType> ii2 = itk::ImageRegionConstIterator<InputImageType>( input2, outputRegionForThread );
	itk::ImageRegionIterator<OutputImageType>     oi  = itk::ImageRegionIterator<OutputImageType>(     output, outputRegionForThread );
	// Scale factor:
	double corr = (double)(m_NumberOfPixels)/(double)(m_NumberOfPixels-1);
	// Iterate through the whole image:
	for( ii1.GoToBegin(), ii2.GoToBegin(), oi.GoToBegin(); !oi.IsAtEnd(); ++ii1, ++ii2, ++oi ){
		double mean   = static_cast<double>( ii1.Get() );
		mean          = mean*mean;
		double sqmean = static_cast<double>( ii2.Get() );
		double cov    = corr*( sqmean - mean );
		cov           = ( cov + 0.1f )/( mean + 0.1f );
		oi.Set(   static_cast<OutputPixelType>( cov )   );		
	}
}

template <class TInputImage, class TOutput>
void ComputeCoefficientOfVariation<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "NumberOfPixels: " << m_NumberOfPixels << std::endl;
}


} // end namespace itk

#endif
