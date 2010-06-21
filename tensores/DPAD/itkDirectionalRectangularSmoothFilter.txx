/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDirectionalRectangularSmoothFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkDirectionalRectangularSmoothFilter_txx
#define _itkDirectionalRectangularSmoothFilter_txx
#include "itkDirectionalRectangularSmoothFilter.h"

#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
DirectionalRectangularSmoothFilter<TInputImage, TOutputImage>
::DirectionalRectangularSmoothFilter()
{
	m_Radius    = 3;
	m_Direction = 0;
}

template <class TInputImage, class TOutputImage>
void DirectionalRectangularSmoothFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	
	// Get pointers to the input and output
	typename Superclass::InputImagePointer inputPtr   = const_cast< TInputImage * >(   this->GetInput()   );
	
	if ( !inputPtr ){ return; }

	if( m_Direction<0 || m_Direction>=TInputImage::ImageDimension )
		itkExceptionMacro( << "The filtering diretion is not appropriate" );
	
	// Get a copy of the input requested region:
	InputImageRegionType inputRequestedRegion = inputPtr->GetRequestedRegion();
	
	// Pad the input requested region by the operator radius
	InputSizeType radius;
	radius.Fill( 0 );
	radius[m_Direction] = m_Radius;
	inputRequestedRegion.PadByRadius( radius );
	// Crop the region to its maximum possible size:
	inputRequestedRegion.Crop(   inputPtr->GetLargestPossibleRegion()   );
	// Set the requested region:
	inputPtr->SetRequestedRegion( inputRequestedRegion );
	return;
}


template< class TInputImage, class TOutputImage>
void DirectionalRectangularSmoothFilter< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Allocate input and output
	OutputImagePointer       output  =  this->GetOutput();
	InputImageConstPointer   input   =  this->GetInput();
	// Compute the corresponding region on the input:
	InputImageRegionType inputRequestedRegion;
	InputSizeType        inputRequestedSize;
	InputIndexType       inputRequestedIndex;
	unsigned long        header = 0;
	unsigned long        tail   = 0;
	for( unsigned int k=0; k<TInputImage::ImageDimension; k++ ){
		if( k==m_Direction ){
			inputRequestedIndex[k] = outputRegionForThread.GetIndex()[k] - m_Radius;
			inputRequestedSize[k]  = outputRegionForThread.GetSize()[k] + 2*m_Radius;
			if( inputRequestedIndex[k]<0 ){  // Initial padding of zeros
				header                 = -inputRequestedIndex[k];
				inputRequestedSize[k] -= header;
				inputRequestedIndex[k] = 0;
			}
			if( inputRequestedIndex[k]+inputRequestedSize[k]>input->GetLargestPossibleRegion().GetSize()[k] ){ // Final padding
				tail                   = inputRequestedIndex[k] + inputRequestedSize[k] - input->GetLargestPossibleRegion().GetSize()[k];
				inputRequestedSize[k] -= tail;
			}
		}
		else{
			inputRequestedIndex[k] = outputRegionForThread.GetIndex()[k];
			inputRequestedSize[k]  = outputRegionForThread.GetSize()[k];
		}
	}
	inputRequestedRegion.SetIndex(   inputRequestedIndex   );
	inputRequestedRegion.SetSize(    inputRequestedSize    );
	// Iterators:
	itk::ImageLinearConstIteratorWithIndex<InputImageType> ii( input,  inputRequestedRegion  ); // Input image
	itk::ImageLinearIteratorWithIndex<OutputImageType>     oi( output, outputRegionForThread ); // Ouput image
	ii.SetDirection( m_Direction );
	oi.SetDirection( m_Direction );
	// Auxiliar parameters:
	double fact    = 1.0/(   2.0*( (double)m_Radius )   +   1.0   );  // Scaling factor
	// Create an auxiliar buffer to store each line:
	unsigned long buffer_size = outputRegionForThread.GetSize()[m_Direction] + 2*m_Radius;
	double* buffer = new double[ buffer_size ];

	for( ii.GoToBegin(), oi.GoToBegin(); !oi.IsAtEnd(); ii.NextLine(), oi.NextLine() ){
		// Auxiliar values:
		//double px   = itk::NumericTraits<double>::Zero;
		double csum = itk::NumericTraits<double>::Zero;
		// Auxiliar counter:
		unsigned long init = 0;
		unsigned long end  = 0;
		// Header padding of zeros:
		for( ; init<header; ++init )
			buffer[init] = itk::NumericTraits<double>::Zero;
		// Actual values from the input image:
		for( ; init<buffer_size-tail; ++init, ++ii )
			buffer[init] = static_cast< double >(   ii.Get()   );
		// Tail padding of zeros:
		for( ; init<buffer_size; ++init )
			buffer[init] = itk::NumericTraits<double>::Zero;
		unsigned long control = 0;
		for( init=0, end=2*m_Radius+1; !oi.IsAtEndOfLine(); ++init, ++end, ++oi ){
			if( control%(2*m_Radius+1) ){
				csum = csum - buffer[init];
				csum = csum + buffer[end];
			}
			else{ // We need to periodically reset the error
				csum = itk::NumericTraits<double>::Zero;
				for( unsigned int k=init; k<end; ++k )
					csum += buffer[k];
			}
			oi.Set(   static_cast<OutputPixelType>( csum*fact )   );
		}
	}

	// Delete the auxiliar buffer:
	delete[] buffer;
}


template <class TInputImage, class TOutputImage>
void DirectionalRectangularSmoothFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Radius: "    << m_Radius    << std::endl;
	os << indent << "Direction: " << m_Direction << std::endl;
}

} // end namespace itk

#endif
