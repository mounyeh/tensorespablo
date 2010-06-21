/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSmoothFieldFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkSmoothFieldFilter_txx
#define _itkSmoothFieldFilter_txx
#include "itkSmoothFieldFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "math.h"

namespace itk
{

template < class TScalar, unsigned int NDimension >
SmoothFieldFilter< TScalar, NDimension >
::SmoothFieldFilter()
{
	m_Sigma.Fill( 2 );
}


template < class TScalar, unsigned int NDimension >
void SmoothFieldFilter< TScalar, NDimension >
::GenerateData( void )
{
	// As well as a smoothing filter:
	SmoothingFilterPointer smoothers[NDimension];
	// Create the first filter in the mini-pipeline:
	smoothers[0] = SmoothingFilterType::New();
	// Set the first filter:
	smoothers[0]->SetRadius( (unsigned int)(m_Sigma[0]) );
	smoothers[0]->SetDirection( 0 );
	smoothers[0]->SetInput( this->GetInput() );
	// Creation of the mini-pipeline
	for( unsigned int d=1; d<NDimension; d++ ){
		smoothers[d] = SmoothingFilterType::New();
		smoothers[d]->SetRadius( (unsigned int)(m_Sigma[d]) );
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
void SmoothFieldFilter< TScalar, NDimension >::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "m_Sigma: " << m_Sigma << std::endl;
}

} // end namespace itk

#endif
