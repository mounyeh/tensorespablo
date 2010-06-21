/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeLCCFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkComputeLCCFilter_txx
#define _itkComputeLCCFilter_txx
#include "itkComputeLCCFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "math.h"

namespace itk
{

template <const unsigned int NDimension>
ComputeLCCFilter<NDimension>
::ComputeLCCFilter()
{
	// this filter requires four input images:
	this->SetNumberOfRequiredInputs(  7 );
	this->SetNumberOfRequiredOutputs( 1 );

	m_Tol = 1e-9;
	m_LCC = 0.0;
}

template <const unsigned int NDimension>
void ComputeLCCFilter<NDimension>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	InputImagePointer  input1  = const_cast< InputImageType* >( this->GetInput1()  );
	InputImagePointer  input2  = const_cast< InputImageType* >( this->GetInput2()  );
	InputImagePointer  input3  = const_cast< InputImageType* >( this->GetInput3()  );
	InputImagePointer  input4  = const_cast< InputImageType* >( this->GetInput4()  );
	InputImagePointer  input5  = const_cast< InputImageType* >( this->GetInput5()  );
	InputImagePointer  input6  = const_cast< InputImageType* >( this->GetInput6()  );
	InputImagePointer  input7  = const_cast< InputImageType* >( this->GetInput7()  );
	OutputImagePointer  output = this->GetOutput();
	
	if (  (!input1) || (!input2) || (!input3) || (!input4) || (!input5)  || (!input6) || (!input7) || (!output)  )
		return;

	// Get a copy of the input requested region (should equal the output requested region)
	InputRegionType inputRequestedRegion = input1->GetRequestedRegion();

	// Set the requested regions for both outputs:
	input1->SetRequestedRegion( inputRequestedRegion );
	input2->SetRequestedRegion( inputRequestedRegion );
	input3->SetRequestedRegion( inputRequestedRegion );
	input4->SetRequestedRegion( inputRequestedRegion );
	input5->SetRequestedRegion( inputRequestedRegion );
	input6->SetRequestedRegion( inputRequestedRegion );
	input7->SetRequestedRegion( inputRequestedRegion );
}

template< const unsigned int NDimension >
void ComputeLCCFilter< NDimension>::BeforeThreadedGenerateData()
{
	// Get the number of threads:
	int n = this->GetNumberOfThreads();
	// Set the appropriate values:
	m_PerThreadLCC.SetSize( n );
	m_PerThreadLCC.Fill( 0.0 );
}

template< const unsigned int NDimension >
void ComputeLCCFilter< NDimension>::AfterThreadedGenerateData()
{
	// Initialize the global LCC:
	m_LCC = 0.0;
	// Sum over each thread:
	for( int n=0; n<m_PerThreadLCC.Size(); n++ )
		m_LCC += m_PerThreadLCC[n];
	// Normalize by the number of pixels:
	typename OutputImageType::SizeType region = this->GetOutput()->GetLargestPossibleRegion().GetSize();
	for( unsigned int d=0; d<NDimension; d++ )
		m_LCC /= ((double)(region[d]));
}

template< const unsigned int NDimension >
void ComputeLCCFilter< NDimension>::ThreadedGenerateData( const OutputRegionType& outputRegionForThread, int threadId )
{
	// Iterators:
	ImageRegionConstIterator< InputImageType > iti1;         // Input1
	ImageRegionConstIterator< InputImageType > iti2;         // Input2
	ImageRegionConstIterator< InputImageType > iti3;         // Input3
	ImageRegionConstIterator< InputImageType > iti4;         // Input4
	ImageRegionConstIterator< InputImageType > iti5;         // Input5
	ImageRegionConstIterator< InputImageType > iti6;         // Input6
	ImageRegionConstIterator< InputImageType > iti7;         // Input7
	ImageRegionIterator< OutputImageType >     ito;          // Output
	
	// Get inputs and output:
	typename InputImageType::ConstPointer input1  = this->GetInput1( );
	typename InputImageType::ConstPointer input2  = this->GetInput2( );
	typename InputImageType::ConstPointer input3  = this->GetInput3( );
	typename InputImageType::ConstPointer input4  = this->GetInput4( );
	typename InputImageType::ConstPointer input5  = this->GetInput5( );
	typename InputImageType::ConstPointer input6  = this->GetInput6( );
	typename InputImageType::ConstPointer input7  = this->GetInput7( );
	typename OutputImageType::Pointer     output  = this->GetOutput( );

	// Create iterators:
	iti1  = ImageRegionConstIterator< InputImageType > ( input1,   outputRegionForThread );
	iti2  = ImageRegionConstIterator< InputImageType > ( input2,   outputRegionForThread );
	iti3  = ImageRegionConstIterator< InputImageType > ( input3,   outputRegionForThread );
	iti4  = ImageRegionConstIterator< InputImageType > ( input4,   outputRegionForThread );
	iti5  = ImageRegionConstIterator< InputImageType > ( input5,   outputRegionForThread );
	iti6  = ImageRegionConstIterator< InputImageType > ( input6,   outputRegionForThread );
	iti7  = ImageRegionConstIterator< InputImageType > ( input7,   outputRegionForThread );
	ito   = ImageRegionIterator< OutputImageType > (     output,   outputRegionForThread );

	// Initialize iterators:
	iti1.GoToBegin( );
	iti2.GoToBegin( );
	iti3.GoToBegin( );
	iti4.GoToBegin( );
	iti5.GoToBegin( );
	iti6.GoToBegin( );
	iti7.GoToBegin( );
	ito.GoToBegin(  );

	// Pixels to evaluate inputs and calculate output:
	InputPixelType  ip1;
	InputPixelType  ip2;
	InputPixelType  ip3;
	InputPixelType  ip4;
	InputPixelType  ip5;
	InputPixelType  ip6;
	InputPixelType  ip7;
	OutputPixelType op;

	// Auxiliar values:
	double sigI;    // Variance of fixed  image
	double sigJ;    // Variance of moving image
	double innerp;  // Inner product of I and J

	// Calculate output:
	while( ! ito.IsAtEnd() ){
		ip1 = iti1.Get(); ip2 = iti2.Get(); ip3 = iti3.Get( ); ip4 = iti4.Get( );
		ip5 = iti5.Get(); ip6 = iti6.Get(); ip7 = iti7.Get( );
		//-------------------------------------------------------------------------
		op[0]       = 0.0;
		op[1]       = 0.0;
		
		sigJ        = ip6 - ip4*ip4;
		sigI        = ip5 - ip3*ip3;
		
		if( (sigI>m_Tol) && (sigJ>m_Tol) ){
			innerp  = ip7 - ip3*ip4;
			op[1]  -= (  innerp/sigJ  )*(  ip2 - ip4  );
			op[1]  += ip1 - ip3;
			sigI   = sqrt( sigI );
			sigJ   = sqrt( sigJ );
			op[1] /= ( sigI*sigJ );
			op[0]  = innerp / ( sigI*sigJ );
		}
		//-------------------------------------------------------------------------
		ito.Set( op );
		++iti1; ++iti2; ++iti3; ++iti4;
		++iti5; ++iti6; ++iti7;
		++ito;
	}
}

/**
 * Standard "PrintSelf" method
 */
template < const unsigned int NDimension >
void ComputeLCCFilter< NDimension >::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "m_Tol: " << m_Tol << std::endl;
  os << indent << "m_LCC: " << m_LCC << std::endl;
}

} // end namespace itk

#endif
