/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDPADNDThomasFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkDPADNDThomasFilter_txx
#define _itkDPADNDThomasFilter_txx

#include "itkDPADNDThomasFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
DPADNDThomasFilter<TInputImage, TOutputImage>::DPADNDThomasFilter()
{
	m_TimeStep      = 2.0f;
	m_DiffusionTerm = static_cast< DiffusionImagePointer >( NULL );
}

template < class TInputImage, class TOutputImage >
void DPADNDThomasFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	// Get pointers to the input and output
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );
	// Set the final requested region
    inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutputImage>
void DPADNDThomasFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
	TOutputImage *out = dynamic_cast<TOutputImage*>(output);
	if ( out )
		out->SetRequestedRegion( out->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutput>
void DPADNDThomasFilter<TInputImage, TOutput>
::GenerateData( )
{
	if( TInputImage::ImageDimension<2 )
		itkExceptionMacro( << "This filter supports only dimenssions above 2" );
	
	// Filters definition:
	ThomasPointer thomas[TInputImage::ImageDimension];
	for( unsigned int k=0; k<TInputImage::ImageDimension; k++ )
		thomas[k] = ThomasType::New();
	AddPointer adder[TInputImage::ImageDimension-1];
	for( unsigned int k=0; k<TInputImage::ImageDimension-1; k++ )
		adder[k] = AddType::New();
	
	// Filters settings:
	for( unsigned int k=0; k<TInputImage::ImageDimension; k++ ){
		thomas[k]->SetTimeStep(      this->GetTimeStep()      );
		thomas[k]->SetDiffusionTerm( this->GetDiffusionTerm() );
		thomas[k]->SetInput(         this->GetInput()         );
		thomas[k]->SetDirection( k );
	}
	// Set the inputs to the adders:
	adder[0]->SetInput1( thomas[0]->GetOutput() );
	adder[0]->SetInput2( thomas[1]->GetOutput() );
	for( unsigned int k=1; k<TInputImage::ImageDimension-1; k++ ){
		adder[k]->SetInput1( adder[k-1]->GetOutput()  );
		adder[k]->SetInput2( thomas[k+1]->GetOutput() );
	}

	// Allocate the output:
	OutputImagePointer outputPtr = this->GetOutput( 0 );
	outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
	outputPtr->Allocate();

	// Graft the output of the filter:
	adder[TInputImage::ImageDimension-2]->GraftOutput( outputPtr );
	adder[TInputImage::ImageDimension-2]->Update();
	adder[TInputImage::ImageDimension-2]->GetOutput()->UpdateOutputInformation();
	// Graft back the output of the filter:
	this->GraftOutput( adder[TInputImage::ImageDimension-2]->GetOutput() );

	return;
}

template <class TInputImage, class TOutput>
void DPADNDThomasFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "TimeStep:      " << m_TimeStep      << std::endl;
	os << indent << "DiffusionTerm: " << m_DiffusionTerm << std::endl;
}

} // end namespace itk

#endif
