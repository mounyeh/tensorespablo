/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkASR1DThomasFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkASR1DThomasFilter_txx
#define _itkASR1DThomasFilter_txx

#include "itkASR1DThomasFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ASR1DThomasFilter<TInputImage, TOutputImage>::ASR1DThomasFilter()
{
	m_TimeStep    = 2.0f;
	m_Alpha       = 2.0f;
	m_TensorField = static_cast< TensorImagePointer >( NULL );
	m_Direction   = 0;
}

template < class TInputImage, class TOutputImage >
void ASR1DThomasFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	// Get pointer to the input
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );	
	if ( !inputPtr ){return;}
	// This filter requires the entire input:
    inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
}

template< class TInputImage, class TOutputImage >
void ASR1DThomasFilter< TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
	TOutputImage *out = dynamic_cast<TOutputImage*>( output );
	if ( out )
		out->SetRequestedRegion( out->GetLargestPossibleRegion() );
}

/** Overwritting of SplitRequestedRegion is needed for the Thomas algorithm to perform appropriately */
//----------------------------------------------------------------------------
template< class TInputImage, class TOutputImage>
int ASR1DThomasFilter< TInputImage, TOutputImage >
::SplitRequestedRegion( int i, int num, OutputImageRegionType& splitRegion )
{
	// Get the output pointer
	OutputImageType * outputPtr = this->GetOutput();
	const typename TOutputImage::SizeType& requestedRegionSize = outputPtr->GetRequestedRegion().GetSize();
	
	int splitAxis;
	typename TOutputImage::IndexType splitIndex;
	typename TOutputImage::SizeType splitSize;
	
	// Initialize the splitRegion to the output requested region
	splitRegion = outputPtr->GetRequestedRegion();
	splitIndex  = splitRegion.GetIndex();
	splitSize   = splitRegion.GetSize();
	
	// Split on the outermost dimension available
	splitAxis = outputPtr->GetImageDimension() - 1;
	while(   ( requestedRegionSize[splitAxis] == 1 ) || (splitAxis == m_Direction)   ){
		--splitAxis;
		if( splitAxis < 0 ){ // cannot split
			itkDebugMacro("  Cannot Split");
			return 1;
		}
	}
	
	// Determine the actual number of pieces that will be generated
	typename TOutputImage::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
	int valuesPerThread = (int)::ceil(range/(double)num);
	int maxThreadIdUsed = (int)::ceil(range/(double)valuesPerThread) - 1;
	
	// Split the region
	if (i < maxThreadIdUsed){
		splitIndex[splitAxis] += i*valuesPerThread;
		splitSize[splitAxis] = valuesPerThread;
	}
	if (i == maxThreadIdUsed){
		splitIndex[splitAxis] += i*valuesPerThread;
		// Last thread needs to process the "rest" dimension being split
		splitSize[splitAxis] = splitSize[splitAxis] - i*valuesPerThread;
	}
	
	// Set the split region ivars
	splitRegion.SetIndex( splitIndex );
	splitRegion.SetSize(  splitSize  );
	
	itkDebugMacro("  Split Piece: " << splitRegion );
	
	return maxThreadIdUsed + 1;
}

template< class TInputImage, class TOutputImage >
void ASR1DThomasFilter< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Iterators typedefs:
	typedef itk::ImageLinearConstIteratorWithIndex< InputImageType  >  InputIteratorType;
	typedef itk::ImageLinearConstIteratorWithIndex< TensorImageType  > TensorIteratorType;
	typedef itk::ImageLinearIteratorWithIndex< OutputImageType >       OutputIteratorType;
	
	// Get pointers to the input and output:
	InputImageConstPointer input  = this->GetInput();
	OutputImagePointer     output = this->GetOutput();
		
	// Test the filtering direction:
	if( (this->GetDirection() >= TInputImage::ImageDimension) || (this->GetDirection()<0) )
		itkExceptionMacro( << "Improper direction selected for filtering" );

	// Size of the region along the filtering direction:
	unsigned int fSize = outputRegionForThread.GetSize()[this->GetDirection()];

	// Buffers for Thomas algorithm:
	double* an = new double[fSize+1];
	double* bn = new double[fSize+1];
	double* yn = new double[fSize+1];
	
	// Iterators creation
	InputIteratorType  inputIterator(    input,           outputRegionForThread   );
	TensorIteratorType tensorIterator(   m_TensorField,   outputRegionForThread   );
	OutputIteratorType outputIterator(   output,          outputRegionForThread   );
	inputIterator.SetDirection(  this->GetDirection() );
	tensorIterator.SetDirection( this->GetDirection() );
	outputIterator.SetDirection( this->GetDirection() );

	// Depending on the filtering direction, get the position of the diagonal element in the difusion tensor:
	unsigned int component = 0;
	for( unsigned int k=0; k<this->GetDirection(); k++ ){
		for( unsigned int d=k; d<TInputImage::ImageDimension; d++ )
			component++;
	}

	// Run through the whole image:
	for( inputIterator.GoToBegin(),tensorIterator.GoToBegin(),outputIterator.GoToBegin(); !outputIterator.IsAtEnd(); inputIterator.NextLine(),tensorIterator.NextLine(),outputIterator.NextLine() ){
		//--------------------------------------------------------------
		// Initialise iterators for the current line:
		inputIterator.GoToBeginOfLine();
		tensorIterator.GoToBeginOfLine();
		// Auxiliar value; it will be useful later on:
		double norm  = 2.0f/(  m_TimeStep*((double)(TInputImage::ImageDimension))  );
		// Auxiliar values for the tensor magnitudes:
		double tprev = m_Alpha; // Constant boundary condition
		double tcurr = tensorIterator.Get()[ component ];
		double tnext = 0.0f;    // Not a valid value
		// In the first iteration, a0=0, b0=1:
		an[0] = 1.0f;
		bn[0] = 0.0f;
		yn[0] = ( inputIterator.Get() )/( (double)(TInputImage::ImageDimension) );
		++tensorIterator; ++inputIterator;
		// The rest of the line:
		unsigned long pos = 1;
		while( !tensorIterator.IsAtEndOfLine() ){
			// Get the corresponding values from the iterators:
			tnext   = tensorIterator.Get()[ component ];
			yn[pos] = ( inputIterator.Get() )/( (double)(TInputImage::ImageDimension) );
			// Recursively get the new alpha (an) and beta (bn) values:
			//     Auxiliar value:
			double divider = ( norm + tprev + 2.0f*tcurr + tnext ) - ( tprev + tcurr )*an[pos-1];
			//     Next alpha and beta:
			an[pos] = (   tcurr + tnext   )/divider;
			bn[pos] = (   norm*(yn[pos-1]) + ( tprev + tcurr )*bn[pos-1]   )/divider;
			// Update the tensor magnitudes:
			tprev = tcurr;
			tcurr = tnext; // tnext is not valid now
			// Increment counters for the next position:
			++tensorIterator;
			++inputIterator;
			++pos;
		}
		// It remains to compute the last position of the an, bn, and yn vectors:
		tnext = m_Alpha; // Constant boundary condition
		// Recursively get the new alpha (an) and beta (bn) values:
		//     Auxiliar value:
		double divider = ( norm + tprev + 2.0f*tcurr + tnext ) - ( tprev + tcurr )*an[pos-1];
		//     Next alpha and beta:
		an[pos] = (   tcurr + tnext   )/divider;
		bn[pos] = (   norm*(yn[pos-1]) + ( tprev + tcurr )*bn[pos-1]   )/divider;
		// Regarding yn[pos+1], we compute it as the last xn value, and we use the same yn vector
		// to store the filtered values:
		yn[pos] = (bn[pos])/( 1.0f - an[pos] );
		// Recursively compute the filtered values:
		for( pos=fSize; pos>0; pos-- )
			yn[pos-1] = an[pos]*yn[pos] + bn[pos];
		// Finally, set the corresponding values to the output iterator
		for( outputIterator.GoToBeginOfLine(); !outputIterator.IsAtEndOfLine(); ++outputIterator, ++pos )
			outputIterator.Set(   static_cast<OutputPixelType>( yn[pos] )   );		
	}
	delete[] an;
	delete[] bn;
	delete[] yn;
}

template <class TInputImage, class TOutput>
void ASR1DThomasFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "TimeStep:    " << m_TimeStep    << std::endl;
	os << indent << "TensorField: " << m_TensorField << std::endl;
}

} // end namespace itk

#endif
