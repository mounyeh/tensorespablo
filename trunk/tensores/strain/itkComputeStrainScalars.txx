/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeStrainScalars.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkComputeStrainScalars_txx
#define _itkComputeStrainScalars_txx
#include "itkComputeStrainScalars.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{



template< class TInputImage, class TOutputImage >
ComputeStrainScalars< TInputImage, TOutputImage >
::ComputeStrainScalars()
{
//	if( TOutputImage::ImageDimension != 3 )
//		itkExceptionMacro( << "This filter is only implemented for 3-D data" );
	m_Scalar = INV;
}



template< class TInputImage, class TOutputImage >
void ComputeStrainScalars< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Allocate images:
	typename OutputImageType::Pointer     output = this->GetOutput();
	typename InputImageType::ConstPointer input  = this->GetInput();
	if( !input || !output )
		itkExceptionMacro( << "Input/Output image has not been set" );
	// Iterators:
	itk::ImageRegionConstIterator<InputImageType> bit = itk::ImageRegionConstIterator<InputImageType>( input, outputRegionForThread );
	itk::ImageRegionIterator<OutputImageType>     it  = itk::ImageRegionIterator<OutputImageType>( output, outputRegionForThread );
	
	EigenValuesArrayType eig;
	double aux1, aux2, aux3;
	unsigned int ord=0;
	switch( m_Scalar ){
		case EIG0:
		case ST0:
			ord = 0;
			break;
		case EIG1:
		case ST1:
			ord = 1;
			break;
		case ST2:
			ord = 2;
			break;
		default:
			ord = 0;
			break;
	}

	// Iterate:
	for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it ){
		OutputPixelType op;
		InputPixelType  ip = bit.Get();
		switch( m_Scalar ){
			case EIG0:
			case EIG1:
				ip.ComputeEigenValues(eig);
				op = static_cast<OutputPixelType>( eig[ord] );
				break;
			case INV:
				op = static_cast<OutputPixelType>( ip.GetInvariant() );
				break;
			case ST0:
			case ST1:
			case ST2:
				op = static_cast<OutputPixelType>( ip[ord] );
				break;
			default:
				op = 0;
				break;
		}
		it.Set( op );
	}
	return;
}


	
} // end namespace itk


#endif
