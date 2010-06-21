/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMaskCCBlockMatchingFilter.txx,v $
  Language:  C++
  Date:      $Date: 2003/12/15 14:13:18 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkMaskCCBlockMatchingFilter_txx
#define _itkMaskCCBlockMatchingFilter_txx

#include "itkMaskCCBlockMatchingFilter.h"

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
MaskCCBlockMatchingFilter<TInputImage1, TInputImage2, TMaskImage>::MaskCCBlockMatchingFilter()
{
	// this filter requires two input images:
	this->SetNumberOfRequiredInputs(  2 );
	this->SetNumberOfRequiredOutputs( 1 );
	
	m_BlockSize.Fill(  5 );
	m_SearchSize.Fill( 5 );
	m_Sample.Fill(     6 );
}


template<class TInputImage1, class TInputImage2, class TMaskImage>
void MaskCCBlockMatchingFilter<TInputImage1, TInputImage2, TMaskImage>::GenerateOutputInformation( void )
{
	// Call the superclass' implementation of this method
	Superclass::GenerateOutputInformation();

	OutputImagePointer       output   =   this->GetOutput();
	InputImage1ConstPointer  input1   =   this->GetInput1();
	InputImage2ConstPointer  input2   =   this->GetInput2();
	
	if (  (!input1) || (!input2) || (!output)  )
		return;

	// Input 1 description
	const typename TInputImage1::SpacingType&   inputSpacing1    = input1->GetSpacing();
	const typename TInputImage1::PointType&     inputOrigin1     = input1->GetOrigin();
	const typename TInputImage1::SizeType&      inputSize1       = input1->GetLargestPossibleRegion().GetSize();
	const typename TInputImage1::IndexType&     inputStartIndex1 = input1->GetLargestPossibleRegion().GetIndex();

	// Output description
	typename       OutputImageType::SpacingType outputSpacing;
	typename       OutputImageType::PointType   outputOrigin;
	typename       OutputImageType::SizeType    outputSize;
	typename       OutputImageType::IndexType   outputStartIndex;

	// Metadata for the output:
	for (unsigned int i = 0; i < TInputImage1::ImageDimension; i++){
		outputSpacing[i]    = inputSpacing1[i] * ( (float)(m_Sample[i]) );
		outputOrigin[i]     = inputOrigin1[i]  + ( inputSpacing1[i] )*( (float)(m_Sample[i]) );
		outputSize[i]       = ( ( (unsigned long)(inputSize1[i]) ) / ( (unsigned long)(m_Sample[i]) ) ) - 1;
		if( ( (unsigned long)(inputSize1[i]) ) % ( (unsigned long)(m_Sample[i]) ) )
			(outputSize[i]) ++;
		if( outputSize[i] < 1 )
			outputSize[i]   = 1;
		outputStartIndex[i] = 0;
	}
	
	// Largest possible region of the output:
	typename OutputImageType::RegionType  outputLargestPossibleRegion;
	outputLargestPossibleRegion.SetSize(  outputSize                  );
	outputLargestPossibleRegion.SetIndex( outputStartIndex            );

	// Set output metadata:
	output->SetLargestPossibleRegion(   outputLargestPossibleRegion );
	output->SetSpacing(                 outputSpacing               );
	output->SetOrigin(                  outputOrigin                );
}





template<class TInputImage1, class TInputImage2, class TMaskImage>
void MaskCCBlockMatchingFilter<TInputImage1, TInputImage2, TMaskImage>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	InputImage1Pointer  input1 = const_cast< TInputImage1* >( this->GetInput1() );
	InputImage2Pointer  input2 = const_cast< TInputImage2* >( this->GetInput2() );
	OutputImagePointer  output = this->GetOutput();
	
	if (  (!input1) || (!input2) || (!output) )
		return;

	// We need to compute the output requested region (size and start index)	
	OutputImageSizeType   outputRequestedRegionSize       = output->GetRequestedRegion().GetSize( );
	OutputImageIndexType  outputRequestedRegionStartIndex = output->GetRequestedRegion().GetIndex();

    SizeType              input1RequestedRegionSize;
	IndexType             input1RequestedRegionStartIndex;

	SizeType2             input2RequestedRegionSize;
	IndexType2            input2RequestedRegionStartIndex;
	
	// And from it, the input requested regions
	for ( unsigned int i = 0; i < TInputImage1::ImageDimension; i++ ){
		input1RequestedRegionSize[i]       =   outputRequestedRegionSize[i]             *        m_Sample[i] ;
		input2RequestedRegionSize[i]       =   outputRequestedRegionSize[i]             *        m_Sample[i] ;
		input1RequestedRegionStartIndex[i] = ( outputRequestedRegionStartIndex[i] + 1 ) * (long)(m_Sample[i]);
		input2RequestedRegionStartIndex[i] = ( outputRequestedRegionStartIndex[i] + 1 ) * (long)(m_Sample[i]);
	}
	
	// Requested region 1 creation:
	RegionType                       input1RequestedRegion            ;
	input1RequestedRegion.SetSize(   input1RequestedRegionSize       );
	input1RequestedRegion.SetIndex(  input1RequestedRegionStartIndex );

	// Requested region 2 creation:
	RegionType2                      input2RequestedRegion            ;
	input2RequestedRegion.SetSize(   input2RequestedRegionSize       );
	input2RequestedRegion.SetIndex(  input2RequestedRegionStartIndex );

	// Pad the input requested regions by the appropriate values:
	input1RequestedRegion.PadByRadius( m_SearchSize + m_BlockSize );
	input2RequestedRegion.PadByRadius( m_BlockSize  );

	// Crop the requested regions by the largest possible regions:
	input1RequestedRegion.Crop( input1->GetLargestPossibleRegion() );
	input2RequestedRegion.Crop( input2->GetLargestPossibleRegion() );

	// Set the requested regions:
	input1->SetRequestedRegion( input1RequestedRegion );
	input2->SetRequestedRegion( input2RequestedRegion );
}


template<class TInputImage1, class TInputImage2, class TMaskImage>
void MaskCCBlockMatchingFilter<TInputImage1, TInputImage2, TMaskImage>
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
	unsigned long p_size = 1;
	unsigned long c_size = 1;	
	for( unsigned long p=0; p<TInputImage1::ImageDimension; p++ ){
		p_size *= (   2*( m_SearchSize[p] ) + 1   );
		c_size *= (   2*( m_BlockSize[p]  ) + 1   );
	}
	// Output pixel:
	OutputPixelType op( p_size );
	/** ######################################################################################################### */
	
	InputImage1ConstPointer input1 = this->GetInput1();
	InputImage2ConstPointer input2 = this->GetInput2();
	OutputImagePointer output      = this->GetOutput();

	IndexType    inputStartIndex1  = input1->GetLargestPossibleRegion().GetIndex();
	IndexType    inputStartIndex2  = input2->GetLargestPossibleRegion().GetIndex();
		
	// Region of the first image that describes the block to match:
	RegionType  region1;
	SizeType    size1;
	IndexType   index1;

	// Region of the second image in wich we must search:
	RegionType2 region2;
	SizeType2   size2;
	IndexType2  index2;

	// Precompute the size of the regions:
	for( unsigned int k=0; k<TInputImage1::ImageDimension; k++ ){
		size1[k] = 2*(m_SearchSize[k]) + 1;
		size2[k] = 1;
	}
	
	// Iterators typedefs:
	typedef itk::ConstNeighborhoodIterator< InputImage1Type >    Input1IteratorType;
	typedef itk::ConstNeighborhoodIterator< InputImage2Type >    Input2IteratorType;
	typedef itk::ImageRegionConstIterator< MaskImageType >       MaskIteratorType;

	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutputIteratorType;

	// Boundary condition:
	ZeroFluxNeumannBoundaryCondition<InputImage2Type> nbc;

	// Iterators creation:
	Input1IteratorType it1;
	Input2IteratorType it2;
	MaskIteratorType   it_mask(this->GetMask(), outputRegionForThread);
	OutputIteratorType it( output, outputRegionForThread );

	unsigned long p_center=(p_size-1)/2;
	
	for ( it.GoToBegin(), it_mask.GoToBegin(); !it.IsAtEnd(); ++it, ++it_mask ){
		if(it_mask.Get()){
			for( p=0; p<p_size; p++ )
				op[p] = 0.0;
			
			op[p_center] = 1.0;
			it.Set(op);
		}else{
			//--------------------------------------------------------------
			// Get the index in the output image:
			typename OutputImageType::IndexType idxo = it.GetIndex();
			//--------------------------------------------------------------
			// Get the respective regions in the input images:
			for( unsigned int k=0; k<TInputImage1::ImageDimension; k++ ){
				index1[k] = (idxo[k] + 1) * (m_Sample[k]) + inputStartIndex1[k] - m_SearchSize[k];
				index2[k] = (idxo[k] + 1) * (m_Sample[k]) + inputStartIndex2[k];
			}
		
			region1.SetSize(  size1   );
			region1.SetIndex( index1 );
			region2.SetSize(  size2   );
			region2.SetIndex( index2 );

			//--------------------------------------------------------------
			// Iterators
			it1 = Input1IteratorType( m_BlockSize, input1, region1 );
			it2 = Input2IteratorType( m_BlockSize, input2, region2 );
		
			it1.OverrideBoundaryCondition( &nbc );
			it2.OverrideBoundaryCondition( &nbc );
		
			unsigned long ptr = 0;
			for( it1.GoToBegin(),it2.GoToBegin(); !it1.IsAtEnd(); ++it1 ){
				// For all possible displacements:
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

				op[ptr] = num / ( sqrt(den1) * sqrt(den2) );			
				++ptr;
			}
		
			// Normalization of the pixel in order to compute prior probabilities:
			mean = 0.0;
			for( p=0; p<p_size; p++ ){
				if( op[p]>0.0 )
					mean += op[p];
				else
					op[p] = 0.0;
			}		
			if( mean>0.0 ){
				for( p=0; p<p_size; p++ ){
					op[p] /= mean;
				}
			}		
			else{
				for( p=0; p<p_size; p++ )
					op[p] = 1.0/((double)p_size);
			}
			it.Set( op );
			//--------------------------------------------------------------
		}
	}
	//---------------------------------------------------------------------------------------------------------------
}


template<class TInputImage1, class TInputImage2, class TMaskImage>
void MaskCCBlockMatchingFilter<TInputImage1, TInputImage2, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	
	os << indent << "m_BlockSize: "   << m_BlockSize  << std::endl;
	os << indent << "m_SearchSize: "  << m_SearchSize << std::endl;
	os << indent << "m_Sample: "      << m_Sample     << std::endl;
}


}// end namespace itk
#endif
