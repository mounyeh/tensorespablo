/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDTIEstimateTensorFilter.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkThresholdToMaskImageFilter_txx
#define _itkThresholdToMaskImageFilter_txx
#include "itkThresholdToMaskImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{



template< class TInputImage, class TOutputImage >
ThresholdToMaskImageFilter< TInputImage, TOutputImage >
::ThresholdToMaskImageFilter()
{
	m_Threshold     = 0.0f;
}



template< class TInputImage, class TOutputImage >
void ThresholdToMaskImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Allocate images:
	typename OutputImageType::Pointer     output = this->GetOutput();
	typename InputImageType::ConstPointer input  = this->GetInput();
	// Iterators:
	itk::ImageRegionConstIterator<InputImageType> bit = itk::ImageRegionConstIterator<InputImageType>( input, outputRegionForThread );
	itk::ImageRegionIterator<OutputImageType>     it  = itk::ImageRegionIterator<OutputImageType>( output, outputRegionForThread );
	// Iterate:
	for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it ){
		OutputPixelType op;
		InputPixelType  ip = bit.Get();
		if(   ip > m_Threshold   )
			op = itk::NumericTraits<OutputPixelType>::One;
		else
			op = itk::NumericTraits<OutputPixelType>::Zero;
		it.Set( op );
	}
	return;
}


	
} // end namespace itk


#endif
