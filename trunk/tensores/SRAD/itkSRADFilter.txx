/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSRADFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkSRADFilter_txx
#define _itkSRADFilter_txx

#include "itkSRADFilter.h"
#include "itkProgressReporter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
SRADFilter<TInputImage, TOutputImage>
::SRADFilter()
{
	m_Iter        = 10;
	m_TimeStep    = 2.0f;
	m_Tensor      = NULL;
	m_ExpConstant = 25.0f;
}

template <class TInputImage, class TOutputImage>
SRADFilter<TInputImage, TOutputImage>
::~SRADFilter()
{
	m_Tensor      = NULL;
}

template < class TInputImage, class TOutputImage >
void SRADFilter<TInputImage, TOutputImage>
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
void SRADFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
	TOutputImage *out = dynamic_cast<TOutputImage*>(output);
	if ( out )
		out->SetRequestedRegion( out->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutputImage>
void SRADFilter<TInputImage, TOutputImage>
::GenerateData( )
{
	// Make sure that the tensor filter exists:
	if(    m_Tensor==static_cast< TensorFilterPointer >( NULL )    ){
		m_Tensor = TensorFilterType::New();
		m_Tensor->SetBeta( this->GetBeta() );
		m_Tensor->SetQ0(   this->GetQ0()   );
	}
	// Allocate the output:
	OutputImagePointer outputPtr = this->GetOutput( 0 );
	outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
	outputPtr->Allocate();
	// Create the other filters:
	InScalePointer  inS    = InScaleType::New();
	OutScalePointer outS   = OutScaleType::New();
	ExpPointer      exp    = ExpType::New();
	LogPointer      log    = LogType::New();
	inS->SetFactor(   1.0f/( this->GetExpConstant() )   );
	inS->SetInput(    this->GetInput()                  );
	exp->SetInput(    inS->GetOutput()                  );
	exp->Update();
	ThomasPointer   thomas = ThomasType::New();
	// Iterate:
	double cumulative_tau = 0.0f;
	for( unsigned int iter=0; iter<m_Iter; iter++ ){
		double q0 = ::exp( -cumulative_tau/6.0f );
		// Set the parameters:
		double current_tau = this->GetTimeStep();
		cumulative_tau += current_tau;
		thomas->SetTimeStep( current_tau );
		// Set the appropriate input...
		if( iter==0 )
			m_Tensor->SetInput( exp->GetOutput()    );
		else
			m_Tensor->SetInput( thomas->GetOutput() );
		// .. and compute the new diffusion tensor
		m_Tensor->SetQ0( q0 );
		m_Tensor->Modified();
		m_Tensor->GetOutput()->SetRequestedRegion(    m_Tensor->GetOutput()->GetLargestPossibleRegion()    );
		m_Tensor->Update();

		// Set the diffusion tensor to the filters
		thomas->SetTensorField(  m_Tensor->GetOutput()  );
		// Set the inputs to the filters:
		if( iter==0 )
			thomas->SetInput( exp->GetOutput()    );
		else
			thomas->SetInput( thomas->GetOutput() );
		// Update the filter and the output information:
		thomas->Modified();
		thomas->Update();
		// Another iteration is done:
		this->InvokeEvent( IterationEvent() );
	}
	log->SetInput(   thomas->GetOutput()    );
	outS->SetFactor( this->GetExpConstant() );
	outS->SetInput(  log->GetOutput()       );
	// Graft the output of the filter:
	outS->GraftOutput( outputPtr );
	outS->Modified();
	outS->Update();
	outS->GetOutput()->UpdateOutputInformation();
	// Graft back the output of the filter:
	this->GraftOutput( outS->GetOutput() );
}

template <class TInputImage, class TOutput>
void SRADFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Iterations:    " << m_Iter     << std::endl;
	os << indent << "Time Step:     " << m_TimeStep << std::endl;
	os << indent << "Tensor filter: " << m_Tensor   << std::endl;
}

} // end namespace itk

#endif
