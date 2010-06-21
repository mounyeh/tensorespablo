/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCoefficientOfVariationToDiffusionTerm.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkCoefficientOfVariationToDiffusionTerm_txx
#define _itkCoefficientOfVariationToDiffusionTerm_txx

#include "itkCoefficientOfVariationToDiffusionTerm.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
CoefficientOfVariationToDiffusionTerm<TInputImage, TOutputImage>
::CoefficientOfVariationToDiffusionTerm()
{
	m_Noise             = 0.0f;
	m_DiffusionTermMode = __Aja__;
}

template< class TInputImage, class TOutputImage>
void CoefficientOfVariationToDiffusionTerm< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Allocate input and output
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input   =  this->GetInput();
	// Iterators:
	itk::ImageRegionConstIterator<InputImageType> bit = itk::ImageRegionConstIterator<InputImageType>( input,  outputRegionForThread );
	itk::ImageRegionIterator<OutputImageType>     it  = itk::ImageRegionIterator<OutputImageType>(     output, outputRegionForThread );
	// Iterate through the whole image:
	for( bit.GoToBegin(), it.GoToBegin(); !bit.IsAtEnd(); ++bit, ++it ){
		double cij = static_cast<double>( bit.Get() );
		// Check how to compute the diffusion term:
		switch( m_DiffusionTermMode ){
		case __Aja__:
			cij        = (   1.0f + ( 1.0f/(cij+1e-9) )   )/(   1.0f + ( 1.0f/(m_Noise+1e-9) )   );
			break;
		case __Yu__:
			cij        = 1.0f   /   (   1.0f   +   (  (cij-m_Noise)/( m_Noise*(1.0f+m_Noise) )  )   );
			break;
		case __Simplified__:
			cij        = ( m_Noise + 1e-9 ) /( cij + 1e-9 );
			break;
		default:
			itkExceptionMacro( << "??? Unknown diffusion term mode" );
			break;
		}
		
		it.Set(   static_cast<OutputPixelType>(cij)   );
	}
}

template <class TInputImage, class TOutput>
void CoefficientOfVariationToDiffusionTerm<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Noise: " << m_Noise << std::endl;
}


} // end namespace itk

#endif
