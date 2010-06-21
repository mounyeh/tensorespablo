/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBayesianRegularizationFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkBayesianRegularizationFilter_txx
#define _itkBayesianRegularizationFilter_txx
#include "itkBayesianRegularizationFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"

#include "math.h"

namespace itk
{

template < const unsigned int NDimension >
BayesianRegularizationFilter< NDimension >
::BayesianRegularizationFilter()
{
  m_Radius.Fill(1);
  m_SearchSize.Fill(0);
  m_Init = false;
  m_Sigma.SetSize( NDimension );
  m_Sigma.Fill( 16.0 );
}

template < const unsigned int NDimension >
void BayesianRegularizationFilter< NDimension >
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
	if( !m_Init )
		Initialize();
    
	// call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();

	// get pointers to the input and output
	typename Superclass::InputImagePointer  inputPtr  = const_cast< InputImageType * >( this->GetInput() );
	typename Superclass::OutputImagePointer outputPtr = this->GetOutput();

	if ( !inputPtr || !outputPtr )
		return;

	// get a copy of the input requested region (should equal the output requested region)
	InputImageRegionType inputRequestedRegion = inputPtr->GetRequestedRegion();

	// pad the input requested region by the operator radius
	inputRequestedRegion.PadByRadius( m_Radius );

	// crop the input requested region at the input's largest possible region
	if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) ){
		inputPtr->SetRequestedRegion( inputRequestedRegion );
		return;
	}
	else{
		// Couldn't crop the region (requested region is outside the largest
		// possible region).  Throw an exception.

		// store what we tried to request (prior to trying to crop)
		inputPtr->SetRequestedRegion( inputRequestedRegion );

		// build an exception
		InvalidRequestedRegionError e(__FILE__, __LINE__);
		OStringStream msg;
		msg << static_cast<const char *>(this->GetNameOfClass()) << "::GenerateInputRequestedRegion()";
		e.SetLocation(msg.str().c_str());
		e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
		e.SetDataObject(inputPtr);
		throw e;
	}
}


template < const unsigned int NDimension >
void BayesianRegularizationFilter< NDimension >::Initialize()
{
	if( m_Sigma.Size()<NDimension )
		itkExceptionMacro(<< "Cannot Initialize BayesianRegularizationFilter before setting the squared variance!");
	// Get the number of possible displacements
	unsigned long w_size = 1;
	for( unsigned int k=0; k<NDimension; k++ )
		w_size *= (2*m_SearchSize[k]+1);
	if( w_size==0 )
		itkExceptionMacro(<< "Cannot Initialize BayesianRegularizationFilter before setting the search size!");
	// Allocate the distances
	m_Distance.SetSize( w_size, w_size );
	// Auxiliar vectors for distance calculation
	double u[NDimension];
	double v[NDimension];
	// Auxiliar values for calculating vector components
	unsigned long num;
	unsigned long den = 1;
	unsigned long den2;
	
	// Calculate the distances
	for( unsigned int p=0; p<NDimension-1; p++ )
		den *= (2*m_SearchSize[p]+1);
	
	for( unsigned long k=0; k<w_size; k++ ){
		// ------------------------------------------------
		// Calculate the coordinates of vector u:
		den2 = den;
		num  = k;
		for( unsigned int p=0; p<NDimension-1; p++ ){
			u[NDimension-p-1] = num/den2;
			num  %= den2;
			den2 /= (2*m_SearchSize[NDimension-p-2]+1);
		}
		u[0] = num;
			
		// ------------------------------------------------
		m_Distance[k][k] = 1.0;
		for( unsigned int l=k+1; l<w_size; l++ ){
			// ------------------------------------------------
			// Calculate the coordinates of vector v:
			den2 = den;
			num  = l;
			for( unsigned int p=0; p<NDimension-1; p++ ){
				v[NDimension-p-1] = num/den2;
				num  %= den2;
				den2 /= (2*m_SearchSize[NDimension-p-2]+1);
			}
			v[0] = num;
			// ------------------------------------------------
			double val = 0.0;
			for( unsigned int p=0; p<NDimension; p++ )
				val += ( u[p] - v[p] ) * ( v[p] - u[p] ) / ( 2.0*(m_Sigma[p]) );
			val = exp( val );
			m_Distance[l][k] = val;
			m_Distance[k][l] = val;
		}
	}
	m_Init = true;
}



template < const unsigned int NDimension >
void BayesianRegularizationFilter< NDimension >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId)
{
	if( !m_Init )
		itkExceptionMacro( << "Filter unInitialized!");

	ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;
	
	// Iterators:
	ConstNeighborhoodIterator<InputImageType> bit;  // Input
	ImageRegionIterator<OutputImageType>      it;   // Output
	
	// Allocate output
	typename OutputImageType::Pointer     output = this->GetOutput();
	typename InputImageType::ConstPointer input  = this->GetInput();
	
	// Find the data-set boundary "faces"
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
	NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> bC;
	faceList = bC(input, outputRegionForThread, m_Radius);
	
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
	
	// support progress methods/callbacks
	ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());  
	
	InputPixelType  ip;
	InputPixelType  ip2;
	OutputPixelType op;
	double vMAX;
	double aux;
	double Norm; // Normalization value

	// Process each of the boundary faces.  These are N-d regions which border
	// the edge of the buffer.
	
	for (fit=faceList.begin(); fit != faceList.end(); ++fit){
		// Iterators:
		bit = ConstNeighborhoodIterator<InputImageType>(m_Radius, input,  *fit );
		it  = ImageRegionIterator<OutputImageType>     (          output, *fit );
		// Boundary Condition
		bit.OverrideBoundaryCondition(&nbc);
		// Begin with the threaded region
		bit.GoToBegin();
		it.GoToBegin();
		
		while ( ! bit.IsAtEnd() ){
			ip = bit.GetCenterPixel();
			op = ip;  // Probability a priori
			// Voxels before central voxel
			for( unsigned long q=0; q<(bit.Size())/2; q++ ){
				ip2 = bit.GetPixel( q );  // take the voxel from the input
				for( unsigned long r=0; r<ip.Size(); r++ ){
					vMAX = 0.0;
					for( unsigned long s=0; s<ip.Size(); s++ ){
						aux = (ip2[s])*(m_Distance[r][s]);
						if( aux>vMAX ) vMAX=aux;
					}
					op[r] *= vMAX;
				}
			}
			// Voxels after central voxel
			for( unsigned long q=(bit.Size())/2+1; q<bit.Size(); q++ ){
				ip2 = bit.GetPixel( q ); // take the voxel from the input
				for( unsigned long r=0; r<ip.Size(); r++ ){
					vMAX = 0.0;
					for( unsigned long s=0; s<ip.Size(); s++ ){
						aux = (ip2[s])*(m_Distance[r][s]);
						if( aux>vMAX ) vMAX=aux;
					}
					op[r] *= vMAX;
				}
			}
			// Normalization of posterior probabilities to sum to 1
			Norm = 0.0;
			for( unsigned int r=0; r<ip.Size(); r++ )
				Norm += op[r];
			for( unsigned int r=0; r<ip.Size(); r++ )
				op[r] /= Norm;
			// Update output pixel
			it.Set( op );
			
			// Update iterators
			++bit;
			++it;
			// Pixel completed
			progress.CompletedPixel();
		}
	}
}


/**
 * Standard "PrintSelf" method
 */
template < const unsigned int NDimension >
void BayesianRegularizationFilter< NDimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Init:       " << m_Init << std::endl;
	os << indent << "SearchSize: " << m_SearchSize << std::endl;
	os << indent << "Variance:   " << m_Sigma << std::endl;
	os << indent << "Distances:  " << m_Distance << std::endl;
}

} // end namespace itk

#endif
