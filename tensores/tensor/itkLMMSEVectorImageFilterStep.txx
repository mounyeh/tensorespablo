/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLMMSEVectorImageFilterStep.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkLMMSEVectorImageFilterStep_txx
#define _itkLMMSEVectorImageFilterStep_txx
#include "itkLMMSEVectorImageFilterStep.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"
#include "vnl/vnl_math.h"

namespace itk
{

/** Constructor */
template <class TInputImage, class TOutputImage>
LMMSEVectorImageFilterStep<TInputImage, TOutputImage>::LMMSEVectorImageFilterStep()
{
	m_Radius.Fill(1);
	m_Channels = 1;
	m_MinimumNumberOfUsedVoxelsFiltering = 1;
	m_UseAbsoluteValue = false;
	m_KeepValue = false;
	m_NoiseVariance = 1.0f;
	m_FirstBaseline = 0;
}

/** The requested input region is larger than the corresponding output, so we need to override this method: */
template <class TInputImage, class TOutputImage>
void LMMSEVectorImageFilterStep<TInputImage, TOutputImage>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
	// Call the superclass' implementation of this method
	Superclass::GenerateInputRequestedRegion();
  
	// Get pointers to the input and output
	InputImagePointer  inputPtr  = const_cast< TInputImage * >( this->GetInput() );
	OutputImagePointer outputPtr = this->GetOutput();
  
	if ( !inputPtr || !outputPtr )
		return;

	// Get a copy of the input requested region (should equal the output
	// requested region)
	InputImageRegionType inputRequestedRegion = inputPtr->GetRequestedRegion();

	// Pad the input requested region by the operator radius
	inputRequestedRegion.PadByRadius( m_Radius );

	// Crop the input requested region at the input's largest possible region
	inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion());
	inputPtr->SetRequestedRegion( inputRequestedRegion );
	return;
}


template< class TInputImage, class TOutputImage>
void LMMSEVectorImageFilterStep< TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Boundary conditions for this filter; Neumann conditions are fine
	ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;	
	// Iterators:
	ConstNeighborhoodIterator<InputImageType> bit;  // Iterator for the input image
	ImageRegionIterator<OutputImageType>      it;   // Iterator for the output image
	// Input and output
	InputImageConstPointer   input   =  this->GetInput();
	OutputImagePointer       output  =  this->GetOutput();
	// Find the data-set boundary "faces"
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType           faceList;
	NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>                                  bC;
	faceList = bC( input, outputRegionForThread, m_Radius );
	typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
	// Auxilair variables to compute statistics for each DWI channel:
	double*       diff                  = new double[m_Channels];
	double*       correction            = new double[m_Channels];
	double*       dSecondAveragedMoment = new double[m_Channels];
	double*       dSquaredMagnitude     = new double[m_Channels];
	unsigned int* iNumberOfUsedVoxels   = new unsigned int[m_Channels];
	double        dFourthAveragedMoment;

	for ( fit=faceList.begin(); fit != faceList.end(); ++fit){ // Iterate through facets
		// Iterators:
		bit = ConstNeighborhoodIterator<InputImageType>(  m_Radius, input, *fit  );
		it  = ImageRegionIterator<OutputImageType>(        output,     *fit      );
		unsigned int neighborhoodSize = bit.Size();
		// Boundary condition:
		bit.OverrideBoundaryCondition( &nbc );
		//===========================================================================================================================
		//===========================================================================================================================
		//===========================================================================================================================
		for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it ){   // Iterate through pixels in the current facet
			//-----------------------------------------------------------------------------------------------------------------------
			// Initiallise the vectors of sample averages for the second order magnitudes, as well as compute the squared value of
			// the central pixel. Initiallise the number of voxels used for the computation of sample statistics in each channel:
			// For the central voxel:
			InputPixelType dMagnitude = bit.GetCenterPixel();
			for ( unsigned int iJ=0; iJ<m_Channels; ++iJ ) { // For each channel
				dSquaredMagnitude[iJ] = dMagnitude[iJ]*dMagnitude[iJ];
				dSecondAveragedMoment[iJ] = 0;
				iNumberOfUsedVoxels[iJ] = 0;
			}
			//-----------------------------------------------------------------------------------------------------------------------
			// With this implementation, we only need to estimate the fourth order moment for the baseline images (or maybe for several
			// baseline images, in the future).
			dFourthAveragedMoment = 0;
			//-----------------------------------------------------------------------------------------------------------------------
			// Compute the sample statistics as the sum oover the neighbourhood:
			for ( unsigned int i = 0; i < neighborhoodSize; ++i ){ // For each voxel in the neighbourhood
				InputPixelType currentPixelValue = bit.GetPixel( i );
				for ( unsigned int iJ=0; iJ<m_Channels; ++iJ ){ // For each channel
					if ( currentPixelValue[iJ]>0 ){ // exactly zero indicates an artifical value filled in by the scanner, maybe make a flag for this test
						iNumberOfUsedVoxels[iJ]++;
						dSecondAveragedMoment[iJ]+= ( currentPixelValue[iJ] * currentPixelValue[iJ] );
					}
				}
				double aux = currentPixelValue[m_FirstBaseline]*currentPixelValue[m_FirstBaseline];
				aux *= aux;
				dFourthAveragedMoment += aux;
			}
			//-----------------------------------------------------------------------------------------------------------------------
			// Compute sample statistics as relative frequencies; if the number of used voxels is too low, keep the central value
			double minSqAvg = itk::NumericTraits<double>::max();
			double maxSqAvg = itk::NumericTraits<double>::min();
			for ( unsigned int iJ=0; iJ<m_Channels; ++iJ ){ // For each DWI channel
				if ( iNumberOfUsedVoxels[iJ]>=m_MinimumNumberOfUsedVoxelsFiltering && dMagnitude[iJ]>0 )
					dSecondAveragedMoment[iJ] /= iNumberOfUsedVoxels[iJ];
				else
					dSecondAveragedMoment[iJ] = dSquaredMagnitude[iJ];
				if( dSecondAveragedMoment[iJ]<minSqAvg )
					minSqAvg = dSecondAveragedMoment[iJ];
				if( dSecondAveragedMoment[iJ]>maxSqAvg )
					maxSqAvg = dSecondAveragedMoment[iJ];
			}
			// The fourth order moment:
			if ( iNumberOfUsedVoxels[m_FirstBaseline]>=m_MinimumNumberOfUsedVoxelsFiltering && dMagnitude[m_FirstBaseline]>0 )
				dFourthAveragedMoment /= iNumberOfUsedVoxels[m_FirstBaseline];
			else
				dFourthAveragedMoment  = dSquaredMagnitude[m_FirstBaseline]*dSquaredMagnitude[m_FirstBaseline];
			//-----------------------------------------------------------------------------------------------------------------------
			// At this point, we have estimates of E{M_i^2,4}, but we need to compute the corresponding estimates for A; however, we
			// compute before the zero-mean input, since dSecondAveragedMoment will be overwriten:
			//   {
			for ( unsigned int iJ=0; iJ<m_Channels; ++iJ )
				diff[iJ] = dSquaredMagnitude[iJ] - dSecondAveragedMoment[iJ];
			//   }
			for ( unsigned int iJ=0; iJ<m_Channels; ++iJ ){
				dSecondAveragedMoment[iJ] -= 2*m_NoiseVariance;
				if( dSecondAveragedMoment[iJ] < 1000*std::numeric_limits<double>::epsilon() )
					dSecondAveragedMoment[iJ] = 1000*std::numeric_limits<double>::epsilon();
			}
			dFourthAveragedMoment -= 8*m_NoiseVariance*( dSecondAveragedMoment[m_FirstBaseline] + m_NoiseVariance );
			if( dFourthAveragedMoment < 1000*std::numeric_limits<double>::epsilon() )
				dFourthAveragedMoment = 1000*std::numeric_limits<double>::epsilon();
			//-----------------------------------------------------------------------------------------------------------------------
			// Now, we have estimates of the moments of A. We have computed as well the difference M - E{M}, that has to be filtered
			// with the inverse of the covariance matrix C_M2M2.
			const unsigned int MAX_ALLOWED_VAR = 1000;
			const float CFACT1 = 5.0f;
			//    -Initiallisation:
			OutputPixelType dFiltered = dMagnitude;
			double normal =  dSecondAveragedMoment[m_FirstBaseline] * dSecondAveragedMoment[m_FirstBaseline];
			//    - Background checking:
			if( normal > 1000*std::numeric_limits<double>::epsilon() ){
				normal        = ( dFourthAveragedMoment - normal ) / normal;
				//    - Variability checking:
				if( normal <= 1000*std::numeric_limits<double>::epsilon() ){     // Variability is extremely slow, so it is likely that an homogeneous region is being filtered
					// In this case, ||C_A2A2|| << ||C_M2M2||, so we simply use the unbiased estimate of the second order moment:
					for( unsigned int iJ=0; iJ<m_Channels; ++iJ )
                        dFiltered[iJ] = dSecondAveragedMoment[iJ];
				}
				else if( normal>MAX_ALLOWED_VAR ){     // Variability is too high, so C_M2M2 is close to singular, and numerical problems may arise; we simply consider tat C_A2M2*(C_M2M2)^(-1) ~= Id
					for( unsigned int iJ=0; iJ<m_Channels; ++iJ )                                                                                                                                                                            
						dFiltered[iJ] = dSecondAveragedMoment[iJ] + diff[iJ];
				}
				else{   // This is the normal case, and should be the one present in the majority of the voxels of the image
					//    - Pre-whitening of the input:
					if( minSqAvg>CFACT1*m_NoiseVariance ) // In this case, the series expansion is convergent, so we may perform the linear correction
						this->CMMInversion( diff, dSecondAveragedMoment, normal, correction, 1 );
					else // The serie expansion is not convergent, and the linear correction is not stable; the aproximation is not accurate, but this corresponds mainly to background pixels, so it is not so important
						this->ComputeInverseMatrix( diff, dSecondAveragedMoment, normal, correction );
					//    - Product with C_A2M2
					//          Scalar product with the vector of second order moments:
					double dp = itk::NumericTraits<double>::Zero;
					for( unsigned int iJ=0; iJ<m_Channels; ++iJ )
						dp += correction[iJ] * dSecondAveragedMoment[iJ];
					normal *= dp;
					//    - Correction of the output value:
					for( unsigned int iJ=0; iJ<m_Channels; ++iJ )
						dFiltered[iJ] = dSecondAveragedMoment[iJ]   +   normal * dSecondAveragedMoment[iJ];
				}
				// Compute the square root of the output, and check if the result is physically consisitent:
				for( unsigned int iJ=0; iJ<m_Channels; ++iJ ){                                                                                                                                                                        
					if( dFiltered[iJ] > 0 )                                                                                                                                                                                       
						dFiltered[iJ] = sqrt( dFiltered[iJ] );                                                                                                                                                              
					else{                                                                                                                                                                                                     
						if ( m_UseAbsoluteValue )                                                                                                                                                                          
							dFiltered[iJ] = sqrt( -dFiltered[iJ] );                                                                                                                                                          
                                                else if( m_KeepValue )                                                                                                                                                                                   
                                                        dFiltered[iJ] = dMagnitude[iJ];                                                                                                                                                                  
                                                else                                                                                                                                                                                                     
                                                        dFiltered[iJ] = 0;                                                                                                                                                                               
                                        }                                                                                                                                                                                                                
                                } 
			}
			else{       // In this case, the second order moment is too small; this is likely to occur in the background
				for( unsigned int iJ=0; iJ<m_Channels; ++iJ ){                                                                                                                                                      
					if ( m_UseAbsoluteValue )                                                                                                                                                                               
						dFiltered[iJ] = sqrt( -dFiltered[iJ] );                                                                                                                                                        
					else if( m_KeepValue )                                                                                                                                                                                   
						dFiltered[iJ] = dMagnitude[iJ];                                                                                                                                                                 
					else                                                                                                                                                                                                    
						dFiltered[iJ] = 0;                                                                                                                                                                              
				}
			}
			//-----------------------------------------------------------------------------------------------------------------------
			// Put the output in place:
			it.Set( dFiltered );
		}
		//===========================================================================================================================
		//===========================================================================================================================
		//===========================================================================================================================
	}
	// Delete previously alloctaed memory:
	delete [] diff;
	delete [] correction;
	delete [] dSecondAveragedMoment;
	delete [] dSquaredMagnitude;
	delete [] iNumberOfUsedVoxels;
}


/** Smart approximate inversion of C_{M^2M^2} (high SNR case)*/
template <class TInputImage, class TOutput>
void LMMSEVectorImageFilterStep<TInputImage, TOutput>
::CMMInversion( const double* measures, const double* squaredAverages, double normal, double* whitened, unsigned int order ) const
{
	// Where:
	//     measures: the squared measurements, which is, the original data (one per channel)
	//     squaredAverages: The vector containing the second order moment for each DWI channel
	//     normal: the variance of the second order moment normalised by the square of the second order moment
	//     whitened: the processed signal, which is, C_MM^(-1)*(M^2-E{M^2})
	//     order: the number of iterations, i.e., the order of Taylor series expansion
	// Auxiliar value to precompute constants:
	normal     = itk::NumericTraits<double>::One / normal; // For convenience
	double aux = 4.0f * m_NoiseVariance * normal;
	// The terms in the inverse matrix:
	double  Ad = aux;
	double* Ai = new double[m_Channels];
	for( unsigned int k=0; k<m_Channels; ++k ){
		Ad   += squaredAverages[k];
		Ai[k] = itk::NumericTraits<double>::One / ( aux * squaredAverages[k] );
	}
	Ad     = -itk::NumericTraits<double>::One / ( aux * Ad );
	// Now, recursively process the output; initiallise w_0 = x
	for( unsigned int k=0; k<m_Channels; ++k )
		whitened[k] = measures[k];
	double cum; // Auxiliar value
	aux *= m_NoiseVariance;
	// Iterate: w_{n+1} = x - D^{-1}w_n
	for( unsigned int o=0; o<order; ++o ){       // If order=0, this loop does nothing!
		// Compute A_d*w
		cum = itk::NumericTraits<double>::Zero;  // Initiallise acumulator
		for( unsigned int k=0; k<m_Channels; ++k )
			cum += whitened[k];
		cum *= Ad;
		// Compute A_i*w
		for( unsigned int k=0; k<m_Channels; ++k )
			whitened[k] = measures[k] - aux*(   Ai[k] * whitened[k] + cum   );
	}
	// Now we have the truncated series of ( Id + D^(-1) )^(-1). It remains to
	// multiplicate by D^(-1):
	// Compute A_d*w
	cum = itk::NumericTraits<double>::Zero; // Initiallise acumulator
	for( unsigned int k=0; k<m_Channels; ++k )
		cum += whitened[k];
	cum *= Ad;
	// Compute A_i*w + A_d*w and normalise
	for( unsigned int k=0; k<m_Channels; ++k )
		whitened[k] = (   Ai[k] * whitened[k] + cum   )*normal;
	// Delete allocated memory:
	delete[] Ai;
	return;
}


/** Smart approximate inversion of C_{M^2M^2} (low SNR extreme)*/
template <class TInputImage, class TOutput>
void LMMSEVectorImageFilterStep<TInputImage, TOutput>
::LowSNRCMMInversion( const double* measures, const double* squaredAverages, double normal, double* whitened, unsigned int order ) const
{
	// Where:
	//     measures: the squared measurements, which is, the original data (one per channel)
	//     squaredAverages: The vector containing the second order moment for each DWI channel
	//     normal: the variance of the second order moment normalised by the square of the second order moment
	//     whitened: the processed signal, which is, C_MM^(-1)*(M^2-E{M^2})
	//     order: the number of iterations, i.e., the order of Taylor series expansion
	// Auxiliar value to precompute constants:
	double aux1 = normal/(4.0f*m_NoiseVariance);
	double aux2 = itk::NumericTraits<double>::One / m_NoiseVariance;
	double* singular = new double[m_Channels]; // Tfe singular part of D
	double* diagonal = new double[m_Channels]; // The diagonal part of D
	// Prepare the matrix:
	for( unsigned int k=0; k<m_Channels; ++k ){
		singular[k] = aux1 * squaredAverages[k];
		diagonal[k] = aux2 * squaredAverages[k];
	}	
	// Now, recursively process the output; initiallise w_0 = x
	for( unsigned int k=0; k<m_Channels; ++k )
		whitened[k] = measures[k];
	// Iterate: w_{n+1} = x - D^{-1}w_n
	for( unsigned int o=0; o<order; ++o ){
		double dp = itk::NumericTraits<double>::Zero; // Auxiliar value to compute the dot product
		for( unsigned int k=0; k<m_Channels; ++k )
			dp += singular[k] * whitened[k];
		for( unsigned int k=0; k<m_Channels; ++k ){
			whitened[k] = measures[k] - (   ( dp + whitened[k] )*diagonal[k]   );
		}
	}
	// Now we have the truncated series of ( Id + D )^(-1). It remains to normalise
	aux2 *= ( 4.0f / m_NoiseVariance );
	for( unsigned int k=0; k<m_Channels; ++k )
		whitened[k] *= aux2;
	// Delete allocated memory:
	delete[] singular;
	delete[] diagonal;
	return;
}

	
/** Matrix inversion; the general case */
template <class TInputImage, class TOutput>
bool LMMSEVectorImageFilterStep<TInputImage, TOutput>
::ComputeInverseMatrix( const double* measures, const double* squaredAverages, double normal, double* whitened ) const
{
	// Compute the matrix to invert
	double** matrix = new double*[m_Channels];
	for( unsigned int j=0; j<m_Channels; ++j )
		matrix[j]    = new double[m_Channels];
	for( unsigned int j=0; j<m_Channels; ++j ){
		matrix[j][j] = normal*squaredAverages[j]*squaredAverages[j] + 4*m_NoiseVariance*(squaredAverages[j] + m_NoiseVariance);
		for( unsigned int k=j+1; k<m_Channels; ++k ){
			matrix[j][k] = normal*squaredAverages[j]*squaredAverages[k];
			matrix[k][j] = matrix[j][k];
		}
	}
	// Compute the independent term:
	double* iterm = new double[m_Channels];
	for( unsigned int j=0; j<m_Channels; ++j )
		iterm[j] = measures[j];
	// For each column col = 1 to m_Channels-1, we need to make zeros in rows from
	// col+1 to m_Channels (note that in C++ array indices are 0-based):
	for( unsigned int col=0; col<m_Channels-1; ++col ){ // For each column
		// We need a non-null element in the position (col,col), in order to
		// accomplish gaussian elimination:
		if( fabs(matrix[col][col]) <= std::numeric_limits<double>::epsilon() ){
			// Bad luck! The element is zero. We need to add a complete row to
			// the row in position c, so that the new element in position (c,c)
			// is not null. Find the first row for which the element (row,col) 
			// is non-zero:
			unsigned int row = col+1;
			while( fabs(matrix[row][col]) <= std::numeric_limits<double>::epsilon() && row<m_Channels ){++row;}
			// If we are not able to find a row satisfying this condition, then
			// the matrix is singular, and this should not be the case; for
			// this reason, we do not perform bound checking, for efficiency. We
			// assume that row is a valid position, and then correct the input
			// and output:
			if( row==m_Channels ){ // Singular matrix!!!
				for( unsigned int j=0; j<m_Channels; ++j ){delete[] matrix[j];}
				delete[] matrix;
				delete[] iterm;
				return false;
			}
			for( unsigned int cc=col; cc<m_Channels; ++cc )
				matrix[col][cc]  += matrix[row][cc];
			iterm[col] += iterm[row];
		}
		// At this point, we have a valid (col,col), element. We scale the whole
		// corresponding col-th row so that the pivoting element is simply 1:
		double scale = itk::NumericTraits<double>::One / matrix[col][col];
		for( unsigned int cc=col; cc<m_Channels; ++cc )
			matrix[col][cc]  *= scale;
		iterm[col] *= scale;
		// Now, we may perform gaussian elimination for each row:
		for( unsigned int row=col+1; row<m_Channels; ++row ){ // For each row
			double scale = matrix[row][col]; // This is the scale, since input[col][col] = 1.
			// Once again, for each column, we add the corresponding scaled
			// version of the pivoting element; however, in the input matrix,
			// values at the left of this column are assumed to be already zero:
			for( unsigned int cc=col; cc<m_Channels; ++cc ) // Only the columns from col
				matrix[row][cc] -= scale * matrix[col][cc];
			iterm[row] -= scale * iterm[col];
			// We have completed this row
		}
		// We have completed this column
	}
	// Now we have an upper-triangular matrix, where all diagonal elements are
	// just 1, except for the last one; Now, we may compute the output in a recursive
	// fashion:
	if( fabs(matrix[m_Channels-1][m_Channels-1]) <= std::numeric_limits<double>::epsilon() ){
		for( unsigned int j=0; j<m_Channels; ++j ){delete[] matrix[j];}
		delete[] matrix;
		delete[] iterm;
		return false;
	}
	whitened[m_Channels-1] = iterm[m_Channels-1] / matrix[m_Channels-1][m_Channels-1]; // The last one
	for( int k=m_Channels-2; k>=0; --k ){ // For each component
		whitened[k] = iterm[k]; // Initiallise
		for( unsigned int j=k+1; j<m_Channels; ++j )
			whitened[k] -= whitened[j] * matrix[k][j];
	}
	// Delete allocated memory:
	for( unsigned int j=0; j<m_Channels; ++j ){delete[] matrix[j];}
	delete[] matrix;
	delete[] iterm;
	// Matrix has been inverted!!
	return true;
}

/** Standard "PrintSelf" method */
template <class TInputImage, class TOutput>
void LMMSEVectorImageFilterStep<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "Radius: "                             << m_Radius                             << std::endl;
	os << indent << "Channels: "                           << m_Channels                           << std::endl;
	os << indent << "UseAbsoluteValue: "                   << m_UseAbsoluteValue                   << std::endl;
	os << indent << "KeepValue: "                          << m_KeepValue                          << std::endl;
	os << indent << "NoiseVariance: "                      << m_NoiseVariance                      << std::endl;
	os << indent << "MinimumNumberOfUsedVoxelsFiltering: " << m_MinimumNumberOfUsedVoxelsFiltering << std::endl;
}

} // end namespace itk

#endif
