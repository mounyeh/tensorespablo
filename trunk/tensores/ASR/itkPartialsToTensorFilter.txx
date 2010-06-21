/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPartialsToTensorFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkPartialsToTensorFilter_txx
#define _itkPartialsToTensorFilter_txx
#include "itkPartialsToTensorFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <math.h>

namespace itk
{

template <class TInputImage, class TOutputImage>
PartialsToTensorFilter<TInputImage, TOutputImage>::PartialsToTensorFilter()
{
	m_Difference = 1.4142f;
	m_Correction = NULL;
}

template< class TInputImage, class TOutputImage>
void PartialsToTensorFilter< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	if( !m_Correction )
		itkExceptionMacro( << "You must set the eigenvalue correction image before updating this filter" );

	// Iterators for each image:
	ImageRegionConstIterator< InputImageType > iit;   // Input
	ImageRegionIterator< OutputImageType >     oit;   // Output
	ImageRegionConstIterator< CorrectionType > cit;   // Correction image

	// Get pointers to the images:
	InputImageConstPointer input   =  this->GetInput();
	OutputImagePointer     output  =  this->GetOutput();
	
	// Create iterators:
	iit = ImageRegionConstIterator< InputImageType >(   input,        outputRegionForThread   );
	oit = ImageRegionIterator< OutputImageType >(       output,       outputRegionForThread   );
	cit = ImageRegionConstIterator< CorrectionType >(   m_Correction, outputRegionForThread   );

	// Initialise:
	iit.GoToBegin();
	oit.GoToBegin();
	cit.GoToBegin();

	// For each pixel in the region:
	while ( ! iit.IsAtEnd() ){
		// The output pixel:
		OutputPixelType op;
		// Compute the principal eigenvalue (the trace of the structure matrix, since all the other
		// eigenvalues must be null)
		double norm = 0.0f;
		double eigenvalue;
		for( unsigned int k=0; k<TInputImage::ImageDimension; k++ )
			norm += (    static_cast<double>( iit.Get()[k] )    )*(    static_cast<double>( iit.Get()[k] )    );
		// Compute the corrected eigenvalue:
		if( norm <= m_Difference )
			eigenvalue = ( 1.0f - (norm*norm)/(m_Difference*m_Difference) );
		else
			eigenvalue = 0.0f;
		// Compute the different components of the diffusion tensor:
		op.Fill( 0.0f );
		if( norm >= 1e-6 ){
			double         tensor;     // Auxiliar value
			unsigned int   pos = 0;    // Position in the sparse matrix that represents the diffusion tensor
			for( unsigned int k=0; k<TInputImage::ImageDimension; k++ ){
				// The diagonal terms have an extra +alpha
				op[pos] = 1.0f;
				for( unsigned int l=k; l<TInputImage::ImageDimension; l++ ){
					tensor   = ( eigenvalue - 1.0f );
					tensor  *= (    static_cast<double>( iit.Get()[k] )    )*(    static_cast<double>( iit.Get()[l] )    );
					tensor  /= norm;
					op[pos] += tensor;
					op[pos] *= (   cit.Get()   );
					pos++;
				}
			}
		}
		else{
			// This is the extreme case of isotropic diffusion:
			unsigned int pos = 0;
			for( unsigned int k=0; k<TInputImage::ImageDimension; k++ ){
				// The diagonal terms are the only not null components of the diffusion tensor, and they all
				// equal alpha
				op[pos] =  cit.Get();
				for( unsigned int l=k; l<TInputImage::ImageDimension; l++ )
					pos++;
			}
		}
		// Set the new pixel value:
		oit.Set( op );
		// Process next pixel:
		++iit; ++oit; ++cit;
	}
}

template <class TInputImage, class TOutput>
void PartialsToTensorFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Difference: " << m_Difference << std::endl;
}



} // end namespace itk

#endif
