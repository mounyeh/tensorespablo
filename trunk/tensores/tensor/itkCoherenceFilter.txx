/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCoherenceFilter.txx,v $
  Language:  C++
  Date:      $Date: 2003/12/15 14:13:18 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkCoherenceFilter_txx
#define _itkCoherenceFilter_txx

#include "itkCoherenceFilter.h"

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
CoherenceFilter<TInputImage1, TInputImage2, TMaskImage>::CoherenceFilter()
{
	// this filter requires two input images:
	this->SetNumberOfRequiredInputs(  2 );
	this->SetNumberOfRequiredOutputs( 1 );
	this->SetNumberOfThreads(1);
	
	m_BlockSize.Fill(  1 );
	m_Label=1;
}


template<class TInputImage1, class TInputImage2, class TMaskImage>
void CoherenceFilter<TInputImage1, TInputImage2, TMaskImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId ){
	// Counter:
	unsigned long p;
	
	unsigned long c_size = 1;	
	for( unsigned long p=0; p<TInputImage1::ImageDimension; p++ ){
		c_size *= (   2*( m_BlockSize[p]  ) + 1   );
	}
	// Output pixel:
	OutputPixelType op;
	
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
	float coh;
	for ( it.GoToBegin(), it1.GoToBegin(), it2.GoToBegin(), it_mask.GoToBegin(); !it.IsAtEnd(); ++it, ++it_mask, ++it1, ++it2 ){
		if(it_mask.Get()!=m_Label){
			op = 0.0;
		}else{
			coh=0;
			for( p=0; p<c_size; p++ ){
				// For all pixels in the block, calculate both mean values:
				coh += fabs(it1.GetCenterPixel() - it1.GetPixel(p));
			}
			op=coh/c_size;
		}
		it.Set( op );
		//--------------------------------------------------------------
	}
	//---------------------------------------------------------------------------------------------------------------
}


template<class TInputImage1, class TInputImage2, class TMaskImage>
void CoherenceFilter<TInputImage1, TInputImage2, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	
	os << indent << "m_BlockSize: "   << m_BlockSize  << std::endl;
	
}


}// end namespace itk
#endif
