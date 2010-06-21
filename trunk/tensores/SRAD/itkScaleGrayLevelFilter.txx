/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScaleGrayLevelFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkScaleGrayLevelFilter_txx
#define _itkScaleGrayLevelFilter_txx

#include "itkScaleGrayLevelFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ScaleGrayLevelFilter<TInputImage, TOutputImage>::ScaleGrayLevelFilter()
{
	m_Factor = 25.0f;
}

template < class TInputImage, class TOutputImage >
void ScaleGrayLevelFilter<TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Allocate input and output
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input   =  this->GetInput();
	itk::ImageRegionIterator<TOutputImage>     bit( output, outputRegionForThread );
	itk::ImageRegionConstIterator<TInputImage> it(  input,  outputRegionForThread );

	for ( bit.GoToBegin(), it.GoToBegin(); !bit.IsAtEnd(); ++bit, ++it){ // Iterate through facets
		double pixel = static_cast<double>(   it.Get()   );
		pixel *= m_Factor;
		bit.Set(   static_cast<OutputPixelType>( pixel )   );
	}
}

template <class TInputImage, class TOutput>
void ScaleGrayLevelFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Factor: " << m_Factor << std::endl;
}

} // end namespace itk

#endif
