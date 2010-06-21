/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkASRFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkASRFilter_txx
#define _itkASRFilter_txx

#include "itkASRFilter.h"
#include "itkProgressReporter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ASRFilter<TInputImage, TOutputImage>::ASRFilter()
{
	m_Iter       = 10;
	m_TimeStep   = 2.0f;
	m_Sigma      = 2.0f;
	m_Difference = 1.4142f;
	m_Beta       = 0.05f;
	m_Tensor     = NULL;
}

template <class TInputImage, class TOutputImage>
ASRFilter<TInputImage, TOutputImage>
::~ASRFilter()
{
	m_Tensor   = NULL;
}

template < class TInputImage, class TOutputImage >
void ASRFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();

	// Get pointer to the input
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );
	if ( !inputPtr )
		return;
	// This filter needs the entire input:
    inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
}

template < class TInputImage, class TOutputImage >
void ASRFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
	TOutputImage *out = dynamic_cast<TOutputImage*>(output);
	if ( out )
		out->SetRequestedRegion( out->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutputImage>
void ASRFilter<TInputImage, TOutputImage>
::GenerateData( )
{
	// Make sure that the tensor filter exists:
	if( m_Tensor==(TensorFilterPointer)NULL ){
		m_Tensor = TensorFilterType::New();
		m_Tensor->SetBeta(       this->GetBeta()       );
		m_Tensor->SetSigma(      this->GetSigma()      );
		m_Tensor->SetDifference( this->GetDifference() );
	}
	// Allocate the output:
	OutputImagePointer outputPtr = this->GetOutput( 0 );
	outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
	outputPtr->Allocate();
	// Create the other filters:
	CRPointer     itoj   = CRType::New();
	ThomasPointer thomas = ThomasType::New();
	// Iterate:
	for( unsigned int iter=0; iter<m_Iter; iter++ ){
		// Set the parameters:
		double current_tau = (  1.0f  -  (double)iter/(double)m_Iter/2.0f  )*( this->GetTimeStep() );
		itoj->SetTimeStep(   current_tau );
		thomas->SetTimeStep( current_tau );
		// Set the appropriate input...
		if( iter==0 )
			m_Tensor->SetInput( this->GetInput()  );
		else
			m_Tensor->SetInput( this->GetOutput() );
		// .. and compute the new diffusion tensor
		m_Tensor->Modified();
		m_Tensor->GetOutput()->SetRequestedRegion(    m_Tensor->GetOutput()->GetLargestPossibleRegion()    );
		m_Tensor->Update();

		// Set the diffusion tensor to the filters
		itoj->SetTensorField(    m_Tensor->GetOutput()  );
		thomas->SetTensorField(  m_Tensor->GetOutput()  );
		// Set the inputs to the filters:
		if( iter==0 )
			itoj->SetInput( this->GetInput() );
		else
			itoj->SetInput( this->GetOutput() );
		thomas->SetInput( itoj->GetOutput() );
		itoj->Modified();
		itoj->Update();
		// Graft the output of the filter:
		thomas->GraftOutput( outputPtr );
		// Update the filter and the output information:
		thomas->Modified();
		thomas->Update();
		thomas->GetOutput()->UpdateOutputInformation();
		// Graft back the output of the filter:
		this->GraftOutput( thomas->GetOutput() );
		// Another iteration is done:
		this->InvokeEvent( IterationEvent() );
	}
}

template <class TInputImage, class TOutput>
void ASRFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Iterations:    " << m_Iter     << std::endl;
	os << indent << "Time Step:     " << m_TimeStep << std::endl;
	os << indent << "Tensor filter: " << m_Tensor   << std::endl;
}

} // end namespace itk

#endif
