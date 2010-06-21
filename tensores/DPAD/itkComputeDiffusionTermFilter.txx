/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeDiffusionTermFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkComputeDiffusionTermFilter_txx
#define _itkComputeDiffusionTermFilter_txx

#include "itkComputeDiffusionTermFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ComputeDiffusionTermFilter<TInputImage, TOutputImage>::ComputeDiffusionTermFilter()
{
	m_Radius.Fill( 2 );
	m_UseMode           = false;
	m_DiffusionTermMode = __Aja__;
}

template < class TInputImage, class TOutputImage >
void ComputeDiffusionTermFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
	// Get pointer to the input
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );
	if ( !inputPtr ){ return; }
	// This filter needs the entire input:
    inputPtr->SetRequestedRegion(   inputPtr->GetLargestPossibleRegion()   );
}

template < class TInputImage, class TOutputImage >
void ComputeDiffusionTermFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
	TOutputImage *out = dynamic_cast<TOutputImage*>(output);
	if ( out )
		out->SetRequestedRegion( out->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutput>
void ComputeDiffusionTermFilter<TInputImage, TOutput>
::GenerateData( )
{
	// Cast the input:
	CastPointer   cast   = CastType::New();
	cast->SetInput(   this->GetInput()   );
	// Compute the square of the image:
	SquarePointer square = SquareType::New();
	square->SetInput(   cast->GetOutput()   );
	// Smooth the input along each direction:
	SmoothPointer smooth1[TInputImage::ImageDimension];
	smooth1[0] = SmoothType::New();
	smooth1[0]->SetDirection(0);
	smooth1[0]->SetRadius( m_Radius[0] );
	smooth1[0]->SetInput( cast->GetOutput() );
	for( unsigned int k=1; k<TInputImage::ImageDimension; ++k ){
		smooth1[k] = SmoothType::New();
		smooth1[k]->SetDirection(k);
		smooth1[k]->SetRadius( m_Radius[k] );
		smooth1[k]->SetInput( smooth1[k-1]->GetOutput() );
	}
	// Smooth the squared input along each direction:
	SmoothPointer smooth2[TInputImage::ImageDimension];
	smooth2[0] = SmoothType::New();
	smooth2[0]->SetDirection(0);
	smooth2[0]->SetRadius( m_Radius[0] );
	smooth2[0]->SetInput( square->GetOutput() );
	for( unsigned int k=1; k<TInputImage::ImageDimension; ++k ){
		smooth2[k] = SmoothType::New();
		smooth2[k]->SetDirection(k);
		smooth2[k]->SetRadius( m_Radius[k] );
		smooth2[k]->SetInput( smooth2[k-1]->GetOutput() );
	}
	// Compute the coefficient of variation:
	unsigned int numberOfPixels = 1;
	for( unsigned int k=0; k<TInputImage::ImageDimension; ++k )
		numberOfPixels *= ( 2*m_Radius[k] + 1 );
	CoefficientOfVariationPointer coeff = CoefficientOfVariationType::New();
	coeff->SetNumberOfPixels( numberOfPixels );
	coeff->SetInput1( smooth1[TInputImage::ImageDimension-1]->GetOutput() );
	coeff->SetInput2( smooth2[TInputImage::ImageDimension-1]->GetOutput() );
	coeff->Update();
	//-----------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------
	// Compute the median value of the coefficient of variation over the entire image:
	//-----------------------------------------------------------------------------------------------
	// Compute the number of elements, and the position of the median:
	unsigned long l=1;
	for( unsigned int k=0; k<TInputImage::ImageDimension; ++k )
		l *= coeff->GetOutput()->GetLargestPossibleRegion().GetSize()[k];
	// Get a pointer to the block of data:
	OutputPixelType* aux = coeff->GetOutput()->GetBufferPointer();
	// Get a copy of the block of data below CV=1
	double*       data     = new double[l]; // The maximum possible value is l
	unsigned long count    = 0;
	for( unsigned long k=0; k<l; ++k ){
		double current_value = static_cast<double>( aux[k] );
		if( current_value<=0.95f ) // Outlier rejection
			data[count++] = current_value;
	}
	double noise;
	if( m_UseMode )
		noise = this->FindMode( data, count );
	else
		noise = this->QuickSort( data, 0, count, count/2, count/2 );
	// Delete the pointer to the data:
	delete[] data;
	//-----------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------
	// And with this number, compute the final value of the diffusion term:
	DiffusionPointer diffusion = DiffusionType::New();
	diffusion->SetInput(   coeff->GetOutput()   );
	diffusion->SetNoise( noise );
	diffusion->SetUseModeDiffusionTerm( this->GetUseModeDiffusionTerm() );
	// Allocate the output:
	OutputImagePointer outputPtr = this->GetOutput( 0 );
	outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
	outputPtr->Allocate();
	// Graft the output of the filter:
	diffusion->GraftOutput( outputPtr );
	// Update the filter and the output information:
	diffusion->Modified();
	diffusion->Update();
	diffusion->GetOutput()->UpdateOutputInformation();
	// Graft back the output of the filter:
	this->GraftOutput( diffusion->GetOutput() );
}

template <class TInputImage, class TOutput>
double ComputeDiffusionTermFilter<TInputImage, TOutput>
::QuickSort( double* list, unsigned long init, unsigned long end, unsigned long initial_divisor, unsigned long target ){
	// Partial quicksort method to partially sort "list" vector until the desired position "target" has been fixed,
	// and then return the list value at this position.

	// Initalization:
	double elem_div = list[initial_divisor];
	long i   = init-1;
	long j   = end;
	int cont = 1;
	if( init>=end )
		return elem_div;
	while( cont ){
		while( list[++i] < elem_div );
		while( list[--j] > elem_div );
		if( i<j ){
			double temp = list[i];
			list[i]     = list[j];
			list[j]     = temp;
		}
		else
			cont = 0;
	}
	double temp           = list[i];
	list[i]               = list[initial_divisor];
	list[initial_divisor] = temp;

	if( i==target )
		return list[i];
	else if( i>target )
		return (   QuickSort( list, init, i-1, i-1, target )   );
	else
		return (   QuickSort( list, i+1,  end, i+1, target )   );
}

template <class TInputImage, class TOutput>
double ComputeDiffusionTermFilter<TInputImage, TOutput>
::FindMode( double* buffer, unsigned long count )
{
	double        mode     = 0.5f;
	unsigned int  R        = 16;
	double        interval = 1.0f;
	// Allocate the buffer for histogram computation:
	const unsigned int bins = 100;
	float histogramB[bins];
	float histogramA[bins];
	while( R>=1 ){
		// Compute the histogram:
		this->ComputeHistogram( buffer, count, 0.0f, mode+interval/2.0f, (float*)histogramB, bins );
		// Create the smoothing filter:
		float* filter = new float[2*R+1];
		this->GaussWin( filter, R );
		// Smooth the histogram and destroy the filter:
		this->SmoothHistogram( (float*) histogramB, (float*) histogramA, bins, filter, 2*R+1 );
		delete[] filter;
		// Find the position of the maximum:
		double maximum   = -1.0f;
		unsigned int pos = bins;
		for( unsigned int k=0; k<bins; ++k ){
			if( histogramA[k]>maximum ){
				maximum = histogramA[k];
				pos     = k;
			}
		}
		// Compute the corresponding value and update the mode:
		double bin_length = ( mode + interval/2.0f )/( (double)bins );
		mode = (   (double)pos + 0.5f   ) * bin_length;
		// Update the width of the smoothing kernel:
		R /= 2;
		// Update the interval to compute the mode:
		interval /= 10.0f;
	}
	return mode;
}

template <class TInputImage, class TOutput>
unsigned long ComputeDiffusionTermFilter<TInputImage, TOutput>
::ComputeHistogram( double* data, unsigned long count, double lower_limit, double upper_limit, float* histogram, unsigned int bins )
{
	// Initialise the histogram with zeros:
	for( unsigned int b=0; b<bins; ++b )
		histogram[b] = 0.0f;
	// The total number of samples inside the desired range:
	unsigned long samples_in_range = 0;
	// The length of each histogram bin:
	double bin_length = ( upper_limit - lower_limit)/((double)bins);
	for( unsigned long k=0; k<count; ++k ){ // For each data sample
		if(   (data[k]>=lower_limit) && (data[k]<=upper_limit)   ){
			// This sample is inside the range:
			samples_in_range++;
			// Compute the position of the sample:
			unsigned int position = (unsigned int)(   ( data[k]-lower_limit )/bin_length   );
			if( position>=bins ){ position=bins-1; }
			// Update the corresponding position in the histogram:
			histogram[position] += 1.0f;
		}
	}
	// Return the total number of samples in range:
	return samples_in_range;
}

template <class TInputImage, class TOutput>
void ComputeDiffusionTermFilter<TInputImage, TOutput>
::SmoothHistogram( float* histogram_before, float* histogram_after, unsigned int bins, float* filter, unsigned int filter_length )
{
	if( bins<filter_length )
		itkExceptionMacro( << "The number of histogram bins must be greater than the filter length" );
	if( 2*(filter_length/2)==filter_length )
		itkExceptionMacro( << "The length of the filter must be an odd number" );
	// Initialise output histogram with zeros:
	for( unsigned int y=0; y<bins; ++y )
		histogram_after[y] = 0.0f;
	// First segment, we need to pad with zeros before
	for( unsigned int y=0; y<filter_length/2; ++y ){
		unsigned int UPLIM = filter_length/2+y+1;
		if( UPLIM>=bins ){ UPLIM=bins-1; }			
		for( unsigned int k=0; k<UPLIM; k++ )
			histogram_after[y] += ( histogram_before[k] )*( filter[filter_length/2-y+k] );
	}
	// Central segment; no zero-padding needed
	for( unsigned int y=filter_length/2; y<bins-filter_length/2; ++y ){
		for( unsigned int k=0; k<filter_length; k++ )
			histogram_after[y] += ( histogram_before[y-filter_length/2+k] )*( filter[k] );
	}
	// Final segment, we need to pad with zeros after
	for( unsigned int y=bins-filter_length/2; y<bins; ++y ){
		unsigned int DOWNLIM = 0;
		if( y-filter_length/2+DOWNLIM<0 ){ DOWNLIM=filter_length/2-y; }
		for( unsigned int k=DOWNLIM; k<filter_length/2+bins-y; k++ )
			histogram_after[y] += ( histogram_before[y-filter_length/2+k] )*( filter[k] );
	}
}

template <class TInputImage, class TOutput>
void ComputeDiffusionTermFilter<TInputImage, TOutput>
::GaussWin( float* filter, unsigned int radius )
{
	// This is a clone of the matlab implementation for gausswin:
	float        alpha = 2.5f;
	unsigned int L     = 2*radius + 1;
	float        N     = (float)L - 1.0f;
	for( unsigned int n=0; n<L; ++n ){
		float arg = 2.0f*alpha*(  (float)n  -  N/2.0f  )/N;
		filter[n] = ::exp(   -arg*arg/2.0f   );
	}
}

template <class TInputImage, class TOutput>
void ComputeDiffusionTermFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Radius: " << m_Radius << std::endl;
}


} // end namespace itk

#endif
