/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDPADFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkDPADFilter_txx
#define _itkDPADFilter_txx

#include "itkDPADFilter.h"
#include "itkProgressReporter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
DPADFilter<TInputImage, TOutputImage>::DPADFilter()
{
	m_Iter              = 10;
	m_TimeStep          = 2.0f;
	m_DiffusionTerm     = NULL;
	m_UseMode           = false;
	m_UseAOS            = true;
	m_DiffusionTermMode = __Aja__;
	m_Radius.Fill( 2 );
}

template <class TInputImage, class TOutputImage>
DPADFilter<TInputImage, TOutputImage>
::~DPADFilter()
{
	m_DiffusionTerm = NULL;
}

template < class TInputImage, class TOutputImage >
void DPADFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	// Get pointer to the input
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );
	if ( !inputPtr ){ return; }
	// This filter needs the entire input:
    inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
}

template < class TInputImage, class TOutputImage >
void DPADFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
	TOutputImage *out = dynamic_cast<TOutputImage*>(output);
	if ( out )
		out->SetRequestedRegion( out->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutputImage>
void DPADFilter<TInputImage, TOutputImage>
::GenerateData( )
{
	// Make sure that the diffusion term filter exists:
	if(  m_DiffusionTerm == (DiffusionTermPointer)NULL  ){
		m_DiffusionTerm = DiffusionTermType::New();
	}
	m_DiffusionTerm->SetRadius( m_Radius );
	m_DiffusionTerm->SetUseMode( m_UseMode );
	m_DiffusionTerm->SetUseModeDiffusionTerm( this->GetUseModeDiffusionTerm() );
	// Create the filters of the pipeline:
	InCastPointer  casti  = InCastType::New();
	ThomasPointer             thomas;
	ForwardDifferencesPointer forwardD;
	if( m_UseAOS ){ thomas   = ThomasType::New(); }
	else{ forwardD = ForwardDifferencesType::New(); }
	OutCastPointer casto  = OutCastType::New();
	// Set the input to the input caster:
	casti->SetInput( this->GetInput() );
	// Iterate:
	for( unsigned int iter=0; iter<m_Iter; iter++ ){
		// Set the parameters:
		double current_tau = (  1.0f  -  (double)iter/(double)m_Iter/2.0f  )*( this->GetTimeStep() );
		if( m_UseAOS ){ thomas->SetTimeStep( current_tau ); }
		else{ forwardD->SetTimeStep( current_tau ); }
		// Set the appropriate input...
		if( iter==0 )
			m_DiffusionTerm->SetInput( casti->GetOutput()  );
		else{
			if( m_UseAOS ){ m_DiffusionTerm->SetInput( thomas->GetOutput() ); }
			else{ m_DiffusionTerm->SetInput( forwardD->GetOutput() ); }
		}
		// .. and compute the new diffusion term
		m_DiffusionTerm->Modified();
		m_DiffusionTerm->GetOutput()->SetRequestedRegion(    m_DiffusionTerm->GetOutput()->GetLargestPossibleRegion()    );
		m_DiffusionTerm->Update();

		// Set the diffusion term to the filters
		if( m_UseAOS ){ thomas->SetDiffusionTerm(  m_DiffusionTerm->GetOutput()  ); }
		else{ forwardD->SetDiffusionTerm(  m_DiffusionTerm->GetOutput()  ); }
		// Set the inputs to the filters:
		if( iter==0 ){
			if( m_UseAOS ){ thomas->SetInput( casti->GetOutput() ); }
			else{ forwardD->SetInput( casti->GetOutput() ); }
		}
		else{
			if( m_UseAOS ){ thomas->SetInput( thomas->GetOutput() ); }
			else{ forwardD->SetInput( forwardD->GetOutput() ); }
		}
		// Set the input to the output caster:
		if( m_UseAOS ){ casto->SetInput( thomas->GetOutput() ); }
		else{ casto->SetInput( forwardD->GetOutput() ); }
		// Allocate the output:
		OutputImagePointer outputPtr = this->GetOutput( 0 );
		outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
		outputPtr->Allocate();
		// Graft the output of the filter:
		casto->GraftOutput( outputPtr );
		// Update the filter and the output information:
		casto->Modified();
		casto->Update();
		casto->GetOutput()->UpdateOutputInformation();
		// Graft back the output of the filter:
		this->GraftOutput( casto->GetOutput() );
		// Another iteration is done:
		this->InvokeEvent( IterationEvent() );
	}
}

template <class TInputImage, class TOutput>
void DPADFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Iterations:     " << m_Iter           << std::endl;
	os << indent << "Time Step:      " << m_TimeStep       << std::endl;
	os << indent << "Diffusion Term: " << m_DiffusionTerm  << std::endl;
}

} // end namespace itk

#endif
