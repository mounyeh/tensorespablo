/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAddVectorsFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkAddVectorsFilter_txx
#define _itkAddVectorsFilter_txx
#include "itkAddVectorsFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "math.h"

namespace itk
{

template < class TScalar, unsigned int NDimension >
AddVectorsFilter< TScalar, NDimension >
::AddVectorsFilter()
{
	// this filter requires two input images:
	this->SetNumberOfRequiredInputs(  2 );
	this->SetNumberOfRequiredOutputs( 1 );
}

template < class TScalar, unsigned int NDimension >
void AddVectorsFilter< TScalar, NDimension >
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	InputImagePointer  input1  = const_cast< InputImageType* >( this->GetInput1() );
	InputImagePointer  input2  = const_cast< InputImageType* >( this->GetInput2() );
	OutputImagePointer  output = this->GetOutput();
	
	if (  (!input1) || (!input2) || (!output) )
		return;

	// Get a copy of the input requested region (should equal the output requested region)
	InputRegionType inputRequestedRegion = input1->GetRequestedRegion();

	// Set the requested regions for both outputs:
	input1->SetRequestedRegion( inputRequestedRegion );
	input2->SetRequestedRegion( inputRequestedRegion );
}

template < class TScalar, unsigned int NDimension >
void AddVectorsFilter< TScalar, NDimension >
::ThreadedGenerateData( const OutputRegionType& outputRegionForThread, int threadId )
{
	// Iterators:
	ImageRegionConstIterator< InputImageType > iti1;        // Input1
	ImageRegionConstIterator< InputImageType > iti2;        // Input2
	ImageRegionIterator< OutputImageType >     ito;         // Output
	
	// Get inputs and output:
	typename InputImageType::ConstPointer  input1 = this->GetInput1();
	typename InputImageType::ConstPointer  input2 = this->GetInput2();
	typename OutputImageType::Pointer      output = this->GetOutput();

	// Create iterators:
	iti1 = ImageRegionConstIterator< InputImageType > ( input1,  outputRegionForThread );
	iti2 = ImageRegionConstIterator< InputImageType > ( input2,  outputRegionForThread );
	ito  = ImageRegionIterator< OutputImageType > (     output,  outputRegionForThread );

	// Initialize iterators:
	iti1.GoToBegin();
	iti2.GoToBegin();
	ito.GoToBegin();

	// Pixels to evaluate inputs and calculate output:
	InputPixelType  ip1;
	InputPixelType  ip2;
	OutputPixelType op;

	// Auxiliar value for the squared norm of the gradient:
	//float norm;

	// Calculate output:
	while( ! ito.IsAtEnd() ){
		ip1 = iti1.Get();
		ip2 = iti2.Get();
		op  = ip1 + ip2;  
		ito.Set( op );
		++iti1;
		++iti2;
		++ito;
	}
}


/**
 * Standard "PrintSelf" method
 */
template < class TScalar, unsigned int NDimension >
void AddVectorsFilter< TScalar, NDimension >::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif
