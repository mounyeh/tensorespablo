/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCorrelationCoefficientFilter.txx,v $
  Language:  C++
  Date:      $Date: 2003/12/15 14:13:18 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkCorrelationCoefficientFilter_txx
#define _itkCorrelationCoefficientFilter_txx

#include "itkCorrelationCoefficientFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"

namespace itk {


template<class TInputImage1, class TInputImage2, class TMaskImage>
CorrelationCoefficientFilter<TInputImage1, TInputImage2, TMaskImage>::CorrelationCoefficientFilter()
{
	// this filter requires two input images:
	this->SetNumberOfRequiredInputs(  2 );
	this->SetNumberOfRequiredOutputs( 1 );
	this->SetNumberOfThreads( 1 );
	m_BlockSize.Fill(  1 );
	m_Label=1;
}


template<class TInputImage1, class TInputImage2, class TMaskImage>
void CorrelationCoefficientFilter<TInputImage1, TInputImage2, TMaskImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId ){
	/** ######################################################################################################### */
	// Counter:
	unsigned long p;
	// Auxiliar values to calculate CC:
	double mu1;
	double mu2;
	double num;
	double den1;
	double den2;
	double mean;	
	// Size of the array of doubles that compose the pixel:
	unsigned long c_size = 1;	
	for( unsigned long p=0; p<TInputImage1::ImageDimension; p++ ){
		c_size *= (   2*( m_BlockSize[p]  ) + 1   );
	}
	// Output pixel:
	OutputPixelType op;
	/** ######################################################################################################### */
	
	InputImage1ConstPointer input1 = this->GetInput1();
	InputImage2ConstPointer input2 = this->GetInput2();
	OutputImagePointer output      = this->GetOutput();
	
	// Iterators typedefs:
	typedef itk::ConstNeighborhoodIterator< InputImage1Type >    Input1IteratorType;
	typedef itk::ConstNeighborhoodIterator< InputImage2Type >    Input2IteratorType;
	typedef itk::ImageRegionConstIterator< MaskImageType >       MaskIteratorType;

	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutputIteratorType;

	// Boundary condition:
	ZeroFluxNeumannBoundaryCondition<InputImage2Type> nbc;

	// Iterators creation:
	Input1IteratorType it1( m_BlockSize, input1, outputRegionForThread );
	Input2IteratorType it2( m_BlockSize, input2, outputRegionForThread );
	MaskIteratorType   it_mask(this->GetMask(), outputRegionForThread);
	OutputIteratorType it( output, outputRegionForThread );

			
	it1.OverrideBoundaryCondition( &nbc );
	it2.OverrideBoundaryCondition( &nbc );

	for ( it.GoToBegin(), it1.GoToBegin(), it2.GoToBegin(), it_mask.GoToBegin(); !it.IsAtEnd(); ++it, ++it_mask, ++it1, ++it2 ){
		if(it_mask.Get()!=m_Label){
			op = 1.0;
			it.Set(op);
		}else{
			mu1 = 0.0;
			mu2 = 0.0;
			
			for( p=0; p<c_size; p++ ){
				// For all pixels in the block, calculate both mean values:
				mu1 += (double)(   it1.GetPixel( p )   ) / (   (double)( c_size )   );
				mu2 += (double)(   it2.GetPixel( p )   ) / (   (double)( c_size )   );
			}
			num  = 0.0;
			den1 = 0.0;
			den2 = 0.0;
			
			for( p=0; p<c_size; p++ ){
				// For all pixels in the block, calculate CC value:
				num  += (   (double)( it1.GetPixel(p) ) - mu1   ) * (   (double)( it2.GetPixel(p) ) - mu2   );
				den1 += (   (double)( it1.GetPixel(p) ) - mu1   ) * (   (double)( it1.GetPixel(p) ) - mu1   );
				den2 += (   (double)( it2.GetPixel(p) ) - mu2   ) * (   (double)( it2.GetPixel(p) ) - mu2   );
			}
			if(den1!=0 && den2!=0){
				op = num / ( sqrt(den1) * sqrt(den2) );
			}else{
				op=1.0;
			}
		}
		
		// Normalization of the pixel in order to compute prior probabilities:
		
		it.Set( op );
		//--------------------------------------------------------------
	}
	//---------------------------------------------------------------------------------------------------------------
}


template<class TInputImage1, class TInputImage2, class TMaskImage>
void CorrelationCoefficientFilter<TInputImage1, TInputImage2, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	
	os << indent << "m_BlockSize: "   << m_BlockSize  << std::endl;
	
}


}// end namespace itk
#endif
