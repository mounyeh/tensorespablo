/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCCBlockMatchingFilter.txx,v $
  Language:  C++
  Date:      $Date: 2003/12/15 14:13:18 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkCCBlockMatchingFilter_txx
#define _itkCCBlockMatchingFilter_txx

#include "itkCCBlockMatchingFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"

namespace itk {


template<class TInputImage1, class TInputImage2>
CCBlockMatchingFilter<TInputImage1, TInputImage2>::CCBlockMatchingFilter()
{
	// this filter requires two input images:
	this->SetNumberOfRequiredInputs(  2 );
	this->SetNumberOfRequiredOutputs( 1 );
	
	m_BlockSize.Fill(  2 );
	m_SearchSize.Fill( 2 );
	m_Sample.Fill(     3 );
	m_Bins     = 64;
	m_MaxValue = 256.0;
	m_Metric   = 7;
}


template<class TInputImage1, class TInputImage2>
void CCBlockMatchingFilter<TInputImage1, TInputImage2>::GenerateOutputInformation( void )
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





template<class TInputImage1, class TInputImage2>
void CCBlockMatchingFilter<TInputImage1, TInputImage2>
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


template<class TInputImage1, class TInputImage2>
void CCBlockMatchingFilter<TInputImage1, TInputImage2>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId ){	
	/** ######################################################################################################### */
	// Counter:
	unsigned long p;
	// Auxiliar values to calculate CR and NMI
	itk::Array< unsigned long >   histFix(  m_Bins );
	itk::Array< unsigned long >   histMov(  m_Bins );
	itk::Array2D< unsigned long > histJoint( m_Bins, m_Bins );
	unsigned int                  posi;
	unsigned int                  posj;
	// Auxiliar entropy values:
	double                        fixedE;
	double                        movingE;
	double                        jointE;
	// Metric values:
	double                        MI;
	double                        NMI;
	// Mean value for normalization purposes:
	double                        mean;
	// Size of the array of doubles that compose the pixel:
	unsigned long p_size = 1;
	unsigned long c_size = 1;	
	for( unsigned long p=0; p<TInputImage1::ImageDimension; p++ ){
		p_size *= (   2*( m_SearchSize[p] ) + 1   );
		c_size *= (   2*( m_BlockSize[p]  ) + 1   );
	}
	// Output pixel:
	OutputPixelType op( p_size );
	itk::Array2D< double > temp( 2, p_size );
	//=============================================================================================================
	double pi;
	double pj;
	double mu;
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
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutputIteratorType;

	// Boundary condition:
	ZeroFluxNeumannBoundaryCondition<InputImage2Type> nbc;

	// Iterators creation:
	Input1IteratorType it1;
	Input2IteratorType it2;
	OutputIteratorType it( output, outputRegionForThread );

	for ( it.GoToBegin(); !it.IsAtEnd(); ++it ){
		// This loop goes trough the whole region that this thread process
		//--------------------------------------------------------------
		// Get the index in the output image:
		typename OutputImageType::IndexType idxo = it.GetIndex();
		//--------------------------------------------------------------
		// Get the respective regions in the input images:
		for( unsigned int k=0; k<TInputImage1::ImageDimension; k++ ){
			index1[k] = (idxo[k] + 1) * (m_Sample[k]) + inputStartIndex1[k] - m_SearchSize[k];
			index2[k] = (idxo[k] + 1) * (m_Sample[k]) + inputStartIndex2[k];
		}

		region1.SetSize(   size1   );
		region1.SetIndex(  index1  );
		region2.SetSize(   size2   );
		region2.SetIndex(  index2  );

		//--------------------------------------------------------------
		// Iterators. Iterator it1 is a neighbourhood iterator whos neighbourhood is the block to match.
		// Its region has size size1, that is, the search size. In other words, it1 goes trough the search
		// vicinity, meanwhile it2, whose size is size2=1x1x...x1 remains fixed. However, we use an iterator
		// in this case to get acess to the pixels in the block to match. Since it1 is defined over input1
		// and it2 over input2, we search, at each subsampled pixel of the moving image (input2), the vector
		// that best represents the assessed displacement, that is, we look for the pixel of the fixed image
		// that best match the current pixel of the moving image, that is, we estimate the inverse transform
		// of that which drives input1 into input2; since the resampling filter we use at last use the inverse
		// transform, we are doing well.
		it1 = Input1IteratorType( m_BlockSize, input1, region1 );
		it2 = Input2IteratorType( m_BlockSize, input2, region2 );
		it1.OverrideBoundaryCondition( &nbc );
		it2.OverrideBoundaryCondition( &nbc );

		//-----------------------------------------------------------------------------------------------------------------------------
		// Precompute statistics of the moving image to reduce computations:
		it2.GoToBegin();
		histMov.Fill( 0 );
		movingE = 0.0;
		for( p=0; p<c_size; p++ ){
			// For every pixel in the moving neighbourhood...
			// ... Get the moving pixel value:
			pj     = (double)(it2.GetPixel(p));
			// ... Get the bin position of the moving pixel:
			posj   = (unsigned int)(      ::floor(   ((double)m_Bins)*pj/m_MaxValue   )      );
			// ... Update histogram count:
			(histMov[posj])++;
		}
		for( unsigned int b=0; b<m_Bins; b++ ){
			if( histMov[b]>0 )                         // The moving entropy
				movingE += (  ::log( (double)(histMov[b]) )  )*( (double)(histMov[b]) );
		}
		movingE = movingE/((double)c_size) - ::log((double)c_size);
		//-----------------------------------------------------------------------------------------------------------------------------
		
		unsigned long ptr = 0;
		for( it1.GoToBegin(); !it1.IsAtEnd(); ++it1 ){
			// This loop goes trough all possible displacements:
			//----------------------------------------------------------------------------------------------------------
			// Initialize histograms and conditioned statistics to compute CR and NMI:
			histFix.Fill(    0  );
			histJoint.Fill(  0  );
			// Initialize entropies:
			fixedE  = 0.0;
			jointE  = 0.0;
			// Initializ metrics:
			MI  = 0.0;
			NMI = 0.0;
			//----------------------------------------------------------------------------------------------------------
						
			for( p=0; p<c_size; p++ ){
				// For every pixel in both neighbourhoods...
				// ... Get the value of the fixed pixel...
				pi           = (double)( it1.GetPixel(p) );
				// ... Get the value of the moving pixel...
				pj           = (double)( it2.GetPixel(p) );
				// ... Get the bin position of the fixed pixel:
				posi         = (unsigned int)(      ::floor(   ((double)m_Bins)*pi/m_MaxValue   )      );
				// ... Get the bin position of the moving pixel:
				posj         = (unsigned int)(      ::floor(   ((double)m_Bins)*pj/m_MaxValue   )      );
				// ... Update histogram counts:
				(histFix[posi])++;
				(histJoint[posj][posi])++;
			}
			for( unsigned int b=0; b<m_Bins; b++ ){
				//-------------------------------------------------------
				// ENTROPY ESTIMATION; to reduce the computations, we do not normalize histogram values, but instead
				// we use histogram counts and normalize entropies later on:		
				if( histFix[b]>0 )                         // The fixed entropy
					fixedE  += (  ::log( (double)(histFix[b]) )  )*( (double)(histFix[b]) );
				for( unsigned int c=0; c<m_Bins; c++ ){    // The joint entropy
					if( histJoint[b][c]>0 )
						jointE += (  ::log( (double)(histJoint[b][c]) )  )*( (double)(histJoint[b][c]) );
				}
			}
			// Entropy normalization:
			fixedE  = fixedE/((double)c_size)  - ::log((double)c_size);
			jointE  = jointE/((double)c_size)  - ::log((double)c_size);
			
			//-----------------------------------------------------------------------------------
			// COMPUTE METRICS
			// [1] Normalized Mutual Information (NMI):
			if( jointE<0.0 )
				NMI = ( fixedE + movingE )/( jointE ) - 1.0;
			else
				NMI = 1.0;
			// [2] Mutual Information (MI). Note that MI equals movingE iff the fixed and the moving
			// block are completelly dependent. Besidies, movingE is constant over all possible displacements, so
			// the normalization is the same for all possible displacements and the computation of probabilities
			// is consistent.
			if( movingE<0.0 )
				MI = ( - fixedE - movingE + jointE ) / ( -movingE );
			//-----------------------------------------------------------------------------------

			//----------------------------------------------------------------------------------------------------------
			temp[0][ptr] = NMI;
			temp[1][ptr] = MI;
			//----------------------------------------------------------------------------------------------------------
			
			++ptr;
		}

		// Normalization of the NMI pixel in order to compute prior probabilities:
		if( m_Metric!=2 ){
			NMI  = 1.0/((double)p_size);
			mean = 0.0;
			for( p=0; p<p_size; p++ ){
				if( temp[0][p]>0.0 )
					mean += temp[0][p];
				else
					temp[0][p] = 0.0;
			}		
			if( mean>0.0 ){
				for( p=0; p<p_size; p++ ){
					temp[0][p] /= mean;
					if( temp[0][p]>NMI )
						NMI = temp[0][p];
				}
			}		
			else{
				for( p=0; p<p_size; p++ )
					temp[0][p] = 1.0/((double)p_size);
			}
		}
		// Normalization of the MI pixel in order to compute prior probabilities:
		if( m_Metric!=1 ){
			MI   = 1.0/((double)p_size);
			mean = 0.0;
			for( p=0; p<p_size; p++ ){
				if( temp[1][p]>0.0 )
					mean += temp[1][p];
				else
					temp[1][p] = 0.0;
			}		
			if( mean>0.0 ){
				for( p=0; p<p_size; p++ ){
					temp[1][p] /= mean;
					if( temp[1][p]>MI )
						MI = temp[1][p];
				}
			}		
			else{
				for( p=0; p<p_size; p++ )
					temp[1][p] = 1.0/((double)p_size);
			}
		}
		
		switch( m_Metric ){
		case 1:
			for( p=0; p<p_size; p++ )
				op[p] = temp[0][p];
			break;
		case 2:
			for( p=0; p<p_size; p++ )
				op[p] = temp[1][p];
			break;
		default:
			//NMI = NMI/(NMI+MI);
			if( NMI>=MI )
				NMI = 1.0;
			MI  = 1.0 - NMI;
			for( p=0; p<p_size; p++ )
				op[p] = NMI*(temp[0][p]) + MI*(temp[1][p]);
			break;
		}
		it.Set( op );
		//--------------------------------------------------------------
	}
	//---------------------------------------------------------------------------------------------------------------
}


template<class TInputImage1, class TInputImage2>
void CCBlockMatchingFilter<TInputImage1, TInputImage2>
::PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	
	os << indent << "m_BlockSize: "   << m_BlockSize  << std::endl;
	os << indent << "m_SearchSize: "  << m_SearchSize << std::endl;
	os << indent << "m_Sample: "      << m_Sample     << std::endl;
}


}// end namespace itk
#endif
