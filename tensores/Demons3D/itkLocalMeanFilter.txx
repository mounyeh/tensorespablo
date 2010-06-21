/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLocalMeanFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkLocalMeanFilter_txx
#define _itkLocalMeanFilter_txx
#include "itkLocalMeanFilter.h"

#include "math.h"

namespace itk
{

template < class TScalar, unsigned int NDimension >
LocalMeanFilter< TScalar, NDimension >
::LocalMeanFilter()
{
	// this filter requires two input images:
	m_Sigma = 4.0;
}


template < class TScalar, unsigned int NDimension >
void LocalMeanFilter< TScalar, NDimension >
::GenerateData( void )
{
	// Get one smoothing filter per image dimension:
	SmoothingFilterPointer smoothers[NDimension];
	// Create the first filter in the mini-pipeline:
	smoothers[0] = SmoothingFilterType::New();
	// Set the first operator:
	smoothers[0]->SetNormalizeAcrossScale( false );
	smoothers[0]->SetSigma( m_Sigma*(this->GetInput()->GetSpacing()[0]) );	
	smoothers[0]->SetZeroOrder();
	smoothers[0]->SetDirection( 0 );
	smoothers[0]->SetInput( this->GetInput() );
	
	// Creation of the mini-pipeline
	for( unsigned int d=1; d<NDimension; d++ ){
		smoothers[d] = SmoothingFilterType::New();
		smoothers[d]->SetNormalizeAcrossScale( false );
		smoothers[d]->SetSigma( m_Sigma*(this->GetInput()->GetSpacing()[d]) );
		smoothers[d]->SetZeroOrder();
		smoothers[d]->SetDirection( d );
		smoothers[d]->SetInput( smoothers[d-1]->GetOutput() );
	}
	
	// Execution of the mini-pipeline:
	smoothers[NDimension-1]->Update();
	
	// Graft the output of the mini-pipeline to this filters' output
	this->GraftOutput( smoothers[NDimension-1]->GetOutput() );
}


/**
 * Standard "PrintSelf" method
 */
template < class TScalar, unsigned int NDimension >
void LocalMeanFilter< TScalar, NDimension >::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "m_Sigma: " << m_Sigma << std::endl;
}

} // end namespace itk

#endif
