/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEigenValueCorrectionFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkEigenValueCorrectionFilter_txx
#define _itkEigenValueCorrectionFilter_txx
#include "itkEigenValueCorrectionFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
EigenValueCorrectionFilter<TInputImage, TOutputImage>::EigenValueCorrectionFilter()
{
	m_Beta    = 0.05f;
	this->SetNumberOfRequiredInputs( 2 );
}

template < class TInputImage, class TOutputImage >
void EigenValueCorrectionFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	this->Superclass::Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	InputImagePointer     input1 =   const_cast< TInputImage * >( this->GetInput1() );
	InputImagePointer     input2 =   const_cast< TInputImage * >( this->GetInput2() );
	OutputImagePointer    output =   this->GetOutput();
	
	if (   !input1   ||   !input2   ||   !output   )
		return;

	// Get a copy of the input requested region (should equal the output requested region)
	input2->SetRequestedRegion(   input1->GetRequestedRegion()   );
    return;
}

template< class TInputImage, class TOutputImage>
void EigenValueCorrectionFilter< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Get pointers to the images:
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input1  =  this->GetInput1();
	InputImageConstPointer   input2  =  this->GetInput2();

	if(  !input1  ||  !input2  ||  !output  )
		itkExceptionMacro( << "Cannot get pointers to all the images" );

	// Iterators for each image:
	ImageRegionConstIterator< InputImageType >     it1( input1, outputRegionForThread );  // Input
	ImageRegionConstIterator< InputImageType >     it2( input2, outputRegionForThread );  // First iteration
	ImageRegionIterator< OutputImageType >         ot(  output, outputRegionForThread );  // Output

	for( it1.GoToBegin(), it2.GoToBegin(), ot.GoToBegin(); !ot.IsAtEnd(); ++it1, ++it2, ++ot ){
		double beta = m_Beta * (   (double)( it1.Get() ) - (double)( it2.Get() )   );
		if( beta<0.0f )
			beta = 1.0f + beta;
		else
			beta = 1.0f - beta;
		if( beta<0.0f ){ beta=0.0f; }
		ot.Set(    static_cast< OutputPixelType >( beta )    );
	}
}

template <class TInputImage, class TOutput>
void EigenValueCorrectionFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Beta: " << m_Beta    << std::endl;
}

} // end namespace itk

#endif
