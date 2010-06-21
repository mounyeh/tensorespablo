/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeCorrectionFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkComputeCorrectionFilter_txx
#define _itkComputeCorrectionFilter_txx
#include "itkComputeCorrectionFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "math.h"

namespace itk
{

template < unsigned int NDimension >
ComputeCorrectionFilter< NDimension >
::ComputeCorrectionFilter()
{
	// this filter requires two input images:
	this->SetNumberOfRequiredInputs(  2 );
	this->SetNumberOfRequiredOutputs( 1 );

	m_Lambda = 1.0;
}

template < unsigned int NDimension >
void ComputeCorrectionFilter< NDimension >
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	InputImage1Pointer  input1 = const_cast< InputImage1Type* >( this->GetInput1() );
	InputImage2Pointer  input2 = const_cast< InputImage2Type* >( this->GetInput2() );
	OutputImagePointer  output = this->GetOutput();
	
	if (  (!input1) || (!input2) || (!output) )
		return;

	// Get a copy of the input requested region (should equal the output requested region)
	Input1RegionType inputRequestedRegion = input1->GetRequestedRegion();

	// Set the requested regions for both outputs:
	input1->SetRequestedRegion( inputRequestedRegion );
	input2->SetRequestedRegion( inputRequestedRegion );
}

template < unsigned int NDimension >
void ComputeCorrectionFilter< NDimension >
::BeforeThreadedGenerateData()
{
	m_Spacing = this->GetInput1()->GetSpacing();
}

template< unsigned int NDimension >
void ComputeCorrectionFilter< NDimension >
::ThreadedGenerateData( const OutputRegionType& outputRegionForThread, int threadId )
{
	// Iterators:
	ImageRegionConstIterator< InputImage1Type > iti1;        // Input1
	ImageRegionConstIterator< InputImage2Type > iti2;        // Input2
	ImageRegionIterator< OutputImageType >      ito;         // Output
	
	// Get inputs and output:
	typename InputImage1Type::ConstPointer input1 = this->GetInput1();
	typename InputImage2Type::ConstPointer input2 = this->GetInput2();
	typename OutputImageType::Pointer      output = this->GetOutput();

	// Create iterators:
	iti1 = ImageRegionConstIterator< InputImage1Type > ( input1,  outputRegionForThread );
	iti2 = ImageRegionConstIterator< InputImage2Type > ( input2,  outputRegionForThread );
	ito  = ImageRegionIterator< OutputImageType > (      output,  outputRegionForThread );

	// Initialize iterators:
	iti1.GoToBegin();
	iti2.GoToBegin();
	ito.GoToBegin();

	// Pixels to evaluate inputs and calculate output:
	Input1PixelType ip1;
	Input2PixelType ip2;
	OutputPixelType op;

	// Auxiliar value for the squared norm of the gradient:
	double norm;

	// Calculate output:
	while( ! ito.IsAtEnd() ){
		ip1   = iti1.Get();
		ip2   = iti2.Get();
		// ip1 is the gradient of the moving image
		// ip2[0] is MATLAB CCaux
		// ip2[1] is MATLAB corr
		// Compute the norm of the gradient of the moving image (note the product by image spacing):
		norm  = 0.0;
		for( unsigned int k=0; k<NDimension; k++ )
			norm += (ip1[k])*(ip1[k])*((double)(m_Spacing[k]))*((double)(m_Spacing[k]));
		// Multiply this value by the square of the correction term to compute the norm of the gradient of the metric:
		norm *= (ip2[1])*(ip2[1]);
		// Add the Levenberg-Marquardt term to this value:
		norm += 4.0*m_Lambda*(ip2[0]);
		// Compute the global correction to the gradient:
		if( norm >1e-2 )
			norm  = 2.0*(ip2[0])/norm;
		else
			norm  = 0.0;
		// In the following loop, we must multiply by the image spacing a couple of times, one to eliminate
		// the dependency of the gradient with the image spacing, one to achieve the warp filter to perform
		// adequately
		for( unsigned int k=0; k<NDimension; k++ )
			op[k] = norm*(ip2[1])*(ip1[k])*((double)(m_Spacing[k]))*((double)(m_Spacing[k]));
		ito.Set( op );
		++iti1;
		++iti2;
		++ito;
	}
}


/**
 * Standard "PrintSelf" method
 */
template < unsigned int NDimension >
void ComputeCorrectionFilter< NDimension >::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "m_Lambda: " << m_Lambda << std::endl;
}

} // end namespace itk

#endif
