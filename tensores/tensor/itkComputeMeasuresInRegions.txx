/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeMeasuresInRegions.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkComputeMeasuresInRegions_txx
#define _itkComputeMeasuresInRegions_txx
#include "itkComputeMeasuresInRegions.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{



template< class TInputImage, class TMaskImage, class TOutputImage >
ComputeMeasuresInRegions< TInputImage, TMaskImage, TOutputImage >
::ComputeMeasuresInRegions()
{
	if( TOutputImage::ImageDimension != 3 )
		itkExceptionMacro( << "This filter is only implemented for 3-D data" );
	m_Scalar = FA;
	m_Label = 0;
}



template< class TInputImage, class TMaskImage, class TOutputImage >
void ComputeMeasuresInRegions< TInputImage, TMaskImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Allocate images:
	typename OutputImageType::Pointer     output = this->GetOutput();
	typename InputImageType::ConstPointer input  = this->GetInput();
	typename MaskImageType::ConstPointer  mask  = this->GetInput2();
	
	if( !input || !output )
		itkExceptionMacro( << "Input/Output image has not been set" );
	// Iterators:
	itk::ImageRegionConstIterator<InputImageType> bit = itk::ImageRegionConstIterator<InputImageType>( input, outputRegionForThread );
	itk::ImageRegionConstIterator<MaskImageType>  mit = itk::ImageRegionConstIterator<MaskImageType>( mask, outputRegionForThread );
	itk::ImageRegionIterator<OutputImageType>     it  = itk::ImageRegionIterator<OutputImageType>( output, outputRegionForThread );
	
	EigenValuesArrayType eig;
	double aux1, aux2, aux3;
	unsigned int ord=0;
	switch( m_Scalar ){
		case EIG0:
		case DC0:
			ord = 0;
			break;
		case EIG1:
		case DC1:
			ord = 1;
			break;
		case EIG2:
		case DC2:
			ord = 2;
			break;
		case DC3:
			ord = 3;
			break;
		case DC4:
			ord = 4;
			break;
		case DC5:
			ord = 5;
			break;
		default:
			ord = 0;
			break;
	}
	// Iterate:
	for( bit.GoToBegin(),it.GoToBegin(), mit.GoToBegin(); !bit.IsAtEnd(); ++bit,++it, ++mit ){
		OutputPixelType op;
		InputPixelType  ip = bit.Get();
		op=NumericTraits<OutputPixelType>::Zero;
		if(mit.Get()==m_Label){
			switch( m_Scalar ){
				case FA:
					op = static_cast<OutputPixelType>( ip.GetFractionalAnisotropy() );
					break;
				case RA:
					op = static_cast<OutputPixelType>( ip.GetRelativeAnisotropy() );
					break;
				case EIG0:
				case EIG1:	
				case EIG2:
					ip.ComputeEigenValues(eig);
					op = static_cast<OutputPixelType>( eig[ord] );
					break;
				case MD:
					op = static_cast<OutputPixelType>( ip.GetMeanDiffusivity() );
					break;
				case DC0:
				case DC1:
				case DC2:
				case DC3:
				case DC4:
				case DC5:
					op = static_cast<OutputPixelType>( ip[ord] );
					break;
				case SC0:
					ip.ComputeShapeCoefficients( aux1, aux2, aux3 );
					op = static_cast<OutputPixelType>( aux1 );
					break;
			case SC1:
				ip.ComputeShapeCoefficients( aux1, aux2, aux3 );
				op = static_cast<OutputPixelType>( aux2 );
				break;
			case SC2:
				ip.ComputeShapeCoefficients( aux1, aux2, aux3 );
				op = static_cast<OutputPixelType>( aux3 );
				break;
			default:
				op = 0;
				break;
			}
		}
		it.Set( op );
	}
	return;
}

	
} // end namespace itk


#endif
