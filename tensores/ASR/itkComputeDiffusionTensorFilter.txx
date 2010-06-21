/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeDiffusionTensorFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkComputeDiffusionTensorFilter_txx
#define _itkComputeDiffusionTensorFilter_txx

#include "itkComputeDiffusionTensorFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ComputeDiffusionTensorFilter<TInputImage, TOutputImage>::ComputeDiffusionTensorFilter()
{
	m_Sigma      = 2.0f;
	m_Difference = 1.4142f;
	m_Beta       = 0.05f;
}

template < class TInputImage, class TOutputImage >
void ComputeDiffusionTensorFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();

	// Get pointer to the input
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );
	if ( !inputPtr )
		return;
	InputImagePointer  firstIt   =   const_cast< TInputImage * >( this->GetFirstIteration() );
	if( !firstIt ){
		this->SetFirstIteration( inputPtr );
		firstIt = inputPtr;
	}
	// This filter needs the entire input:
    inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
	firstIt->SetRequestedRegion(  inputPtr->GetLargestPossibleRegion() );
}

template < class TInputImage, class TOutputImage >
void ComputeDiffusionTensorFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
	TOutputImage *out = dynamic_cast<TOutputImage*>(output);
	if ( out )
		out->SetRequestedRegion( out->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutput>
void ComputeDiffusionTensorFilter<TInputImage, TOutput>
::GenerateData( )
{
	InputImagePointer firstIt = const_cast< TInputImage * >( this->GetFirstIteration() );
	if( !firstIt )
		this->SetFirstIteration( this->GetInput() );

	// The three auxiliar filters:
	GradientFilterPointer gradient = GradientFilterType::New();
	TensorFilterPointer   tensor   = TensorFilterType::New();
	CorrectionPointer     correct  = CorrectionType::New();

	// Set parameters:
	gradient->SetSigma( m_Sigma );
	tensor->SetDifference( m_Difference );
	correct->SetBeta(  this->GetBeta()  );
	
	// Set the inputs to the filters:
	gradient->SetInput( this->GetInput()      );
	tensor->SetInput(   gradient->GetOutput() );
	correct->SetInput1( this->GetFirstIteration() );
	correct->SetInput2( this->GetInput()          );
	correct->Update();
	tensor->SetCorrection( correct->GetOutput() );

	// Allocate the output:
	OutputImagePointer outputPtr = this->GetOutput( 0 );
	outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
	outputPtr->Allocate();

	// Graft the output of the filter:
	tensor->GraftOutput( outputPtr );

	// Update the filter and the output information:
	tensor->Modified();
	tensor->Update();
	tensor->GetOutput()->UpdateOutputInformation();

	// Graft back the output of the filter:
	this->GraftOutput( tensor->GetOutput() );
}

template <class TInputImage, class TOutput>
void ComputeDiffusionTensorFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Sigma:      " << m_Sigma      << std::endl;
	os << indent << "Difference: " << m_Difference << std::endl;
	os << indent << "Beta: "       << m_Beta       << std::endl;
}

} // end namespace itk

#endif
