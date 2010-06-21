/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMAPDisplacementFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkMAPDisplacementFilter_txx
#define _itkMAPDisplacementFilter_txx
#include "itkMAPDisplacementFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"

#include "math.h"

namespace itk
{

template < class OutputPointType, unsigned int NDimension >
MAPDisplacementFilter< OutputPointType, NDimension >::MAPDisplacementFilter()
{
	m_Init = false;
}

template < class OutputPointType, unsigned int NDimension >
void MAPDisplacementFilter< OutputPointType, NDimension >
::Initialize( SpacingType spacing, InputSizeType searchsize )
{
	// Update member variables
	m_Spacing    = spacing;
	m_SearchSize = searchsize;
	m_Init       = true;
	// Get the number of possible displacements
	unsigned long w_size = 1;
	for( unsigned int k=0; k<NDimension; k++ )
		w_size *= (2*m_SearchSize[k]+1);
	
	// Allocate coordinate matrix
	m_Displacement.SetSize(NDimension,w_size);
	
	// Auxiliar values for calculating vector components
	long num;
	long den = 1;
	long den2;
	
	for( unsigned int p=0; p<NDimension-1; p++ )
		den *= (2*m_SearchSize[p]+1);
			
	for( long k=0; k<w_size; k++ ){
		den2 = den;
		num  = k;
		for( long p=0; p<NDimension-1; p++ ){
			m_Displacement[NDimension-p-1][k] = ( (m_Spacing[NDimension-p-1]) * ( (double)( num/den2 - (long)(m_SearchSize[NDimension-p-1]) ) ) );
			num  %= den2;
			den2 /= (2*m_SearchSize[NDimension-p-2]+1);
		}
		m_Displacement[0][k] = ( (m_Spacing[0]) * (double)( num - (long)(m_SearchSize[0]) ) );
	}
}







template < class OutputPointType, unsigned int NDimension >
void MAPDisplacementFilter< OutputPointType, NDimension >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	if( !m_Init )
		itkExceptionMacro( << "Filter uninitialized!");

	// Iterators:
	ImageRegionConstIterator<InputImageType> iti;        // Input
	ImageRegionIterator<OutputImageType>     ito;        // Output
	
	// Allocate output
	typename InputImageType::ConstPointer input  = this->GetInput();
	typename OutputImageType::Pointer     output = this->GetOutput();
	
	// support progress methods/callbacks
	ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());  
	
	InputPixelType  ip;
	OutputPixelType op;

	iti = ImageRegionConstIterator<InputImageType> ( input,  outputRegionForThread );
	ito = ImageRegionIterator<OutputImageType>(      output, outputRegionForThread );

	iti.GoToBegin();
	ito.GoToBegin();

	double        max;
	unsigned long max_idx;

	while( ! iti.IsAtEnd() ){
		ip = iti.Get();

		max     = 0.0;
		max_idx = 0;

		for( unsigned long m=0; m<ip.Size(); m++ ){
			if( ip[m] > max ){
				max     = ip[m];
				max_idx = m;
			}
		}

		for( unsigned int n=0; n<NDimension; n++ )
			op[n] = m_Displacement[n][max_idx];
					
		// Set the value at the output:
		ito.Set( op );
		// Update iterators:
		++iti;
		++ito;
		// Pixel processed:
		progress.CompletedPixel();
	}
}


/**
 * Standard "PrintSelf" method
 */
template < class OutputPointType, unsigned int NDimension >
void MAPDisplacementFilter< OutputPointType, NDimension >::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif
