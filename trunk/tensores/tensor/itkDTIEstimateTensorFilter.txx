/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDTIEstimateTensorFilter.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDTIEstimateTensorFilter_txx
#define _itkDTIEstimateTensorFilter_txx
#include "itkDTIEstimateTensorFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "math.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
DTIEstimateTensorFilter< TInputImage, TOutputImage >
::DTIEstimateTensorFilter()
{
	m_Channels      = 0;
	m_NBaselines    = 0;
	m_B             = -100.0f;
	m_Iterations    = 0;
	m_X             = LSMatrixType( 0 );
	m_PseudoInverse = LSMatrixType( 0 );
	m_Inverse.Fill( 0.0f );
	m_Indicator     = IndicatorType( 0 );
	m_Baselines     = IndicatorType( 0 );
	m_Threshold     = 0.0f;
	m_Mask          = NULL;
	m_Features      = NULL;
	m_T2            = NULL;
	m_ComputeT2     = false;
}

	
template< class TInputImage, class TOutputImage >
void DTIEstimateTensorFilter< TInputImage, TOutputImage >
::Configure()
{
	if(   strcmp( this->GetInput()->GetNameOfClass(), "DWImages" )   )
	   itkExceptionMacro( << "This method may only be used with inputs of type itk::DWImages" );
	//DWIType* dwi = dynamic_cast<DWIType*>(   const_cast<InputImageType*>( this->GetInput() )   );
	InputImageType* dwi = const_cast<InputImageType*>( this->GetInput() );
	this->SetChannels( dwi->GetNumDWImages() );
	this->SetNBaselines( dwi->GetNumBaselines() );
	this->SetDWIChannels( dwi->GetIndexesDWImages() );
	this->SetBaselines( dwi->GetIndexesBaselines() );
	this->SetB( dwi->GetB_Value() );
	typename DWIType::DiffusionDirectionsType directions;
	dwi->GetDiffusionDirections( directions );
	std::cout << " using directions : ";
	for( unsigned int k=0; k<this->GetChannels(); ++k ) {
	    std::cout << directions[k][0] << " "<< directions[k][1] << " " << directions[k][2] << std::endl;
		this->AddGradientDirection( directions[k] );}
	return;
}


	
template< class TInputImage, class TOutputImage >
void DTIEstimateTensorFilter< TInputImage, TOutputImage >
::ComputeThreshold()
{
	StatisticsFilterPointer stats = StatisticsFilterType::New();
	ThresholdFilterPointer  thres = ThresholdFilterType::New();
	stats->SetRadius(2);
	stats->SetChannels( m_NBaselines + m_Channels );
	stats->SetUseNeighborhoodGradients();
	stats->SetIndicator( m_Baselines );
	stats->SetInput( this->GetInput() );
	stats->Update();
	thres->SetMin( stats->GetMin() );
	thres->SetMax( stats->GetMax() );
	thres->SetBins( 2048 );
	thres->SetW( 2.0f );
	thres->SetInput( stats->GetOutput() );
	thres->Update();
	m_Threshold = thres->GetThreshold();
	m_Features = thres->GetOutput();
	stats = NULL;
	thres = NULL;
}

template< class TInputImage, class TOutputImage >
void DTIEstimateTensorFilter< TInputImage, TOutputImage >
::ComputeMask()
{
	if( !m_Features ){
		double th = m_Threshold;
		this->ComputeThreshold();
		if( th>1.0f )
			m_Threshold = th;
	}
	MaskFilterPointer mask = MaskFilterType::New();
	mask->SetInput( m_Features );
	mask->SetThreshold(m_Threshold);
	mask->Update();
	m_Mask = mask->GetOutput();
	m_Features = NULL;
	mask = NULL;
}

template< class TInputImage, class TOutputImage >
void DTIEstimateTensorFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId )
{
	// Allocate images:
	typename OutputImageType::Pointer     output = this->GetOutput();
	typename InputImageType::ConstPointer input  = this->GetInput();
	// Iterators:

	itk::ImageRegionConstIterator<InputImageType> bit = itk::ImageRegionConstIterator<InputImageType>( input, outputRegionForThread );
	itk::ImageRegionConstIterator<MaskImageType>  mit;
	itk::ImageRegionIterator<T2ImageType>         t2;
	itk::ImageRegionIterator<OutputImageType>     it  = itk::ImageRegionIterator<OutputImageType>( output, outputRegionForThread );
	if( m_Mask != (MaskImagePointer)NULL ){
		mit  = itk::ImageRegionConstIterator<MaskImageType>( m_Mask, outputRegionForThread );
		mit.GoToBegin();
	}
	if( m_ComputeT2 ){
		t2  = itk::ImageRegionIterator<T2ImageType>( m_T2, outputRegionForThread );
		t2.GoToBegin();
	}
	// Auxiliar vector to store the input temporarilly:
	double* dwi = new double[m_Channels+1];
	double top[DTI_ESTIMATE_DIM];
	// Auxiliar constants to avoid numeric issues:
	//const double eps  = 0.25f * m_Threshold;
	const double eps  = 1.0f;
	const double epsL = ::log( eps );
	// Baseline:
	double bnorm = 1.0f/(double)m_NBaselines;
	// Iterate:
	for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it ){
		OutputPixelType op;
		// Check the mask:
		if(   !m_Mask   ||   mit.Get()>0 ){
			InputPixelType  ip = bit.Get();
			// Compute the logarithm of the input DWI components:
			for( unsigned int c=0; c<m_Channels; ++c ){
				if( ip[ m_Indicator[c] ] > eps )
					dwi[c] = ::log( ip[ m_Indicator[c] ] );
				else
					dwi[c] = epsL;
			}
			// Compute the logarithm of the [averaged] baseline:
	                double baseline = itk::NumericTraits<double>::Zero;
        	        for( unsigned int b=0; b<m_NBaselines; ++b ){
                	        if( ip[ m_Baselines[b] ] > eps )
                        	        baseline += ::log( ip[ m_Baselines[b] ] );
                       		else
                                	baseline += epsL;
                	}
                	dwi[m_Channels] = baseline*bnorm;
			// Independently on the number of iterations, we must solve the least
			// squares problem; fortunately, it reduces to a simple matrix product
			// with the pseudo-inverse matrix that we have precomputed:	
			for( unsigned int k=0; k<DTI_ESTIMATE_DIM; ++k ){
				top[k] = itk::NumericTraits<double>::Zero;
				for( unsigned int c=0; c<m_Channels+1; ++c )
					top[k] += m_PseudoInverse[c][k] * dwi[c];
			}
			// At this point, we have estimated the components of the tensor with a
			// simple LS approach; however, the LS solution is not the one minimising
			// the mean squared error if the variance of each component is not the same.
			// From the paper by Salvador et al. we know that the variance is
			// approximately proportional to the exponential magnitude exp(-b·g'Dg)
			// when Rician noise is present, so we may use this result to further
			// refine the solution with a WLS approach.
			for( unsigned int iter=0; iter<m_Iterations; ++iter ){
				// Compute weighting factors:
				double  wnorm   = itk::NumericTraits<double>::Zero;
				double* weights = new double[m_Channels+1];
				for( unsigned int c=0; c<m_Channels; ++c ){
					weights[c] = itk::NumericTraits<double>::Zero;
					for( unsigned int k=0; k<DTI_ESTIMATE_DIM; ++k )
						weights[c] += m_X[c][k] * top[k];
					weights[c]  = ::exp( - 2 * weights[c] );
					wnorm      += weights[c];
				}
				weights[m_Channels]  = ::exp( -2.0f * m_X[m_Channels][0] * top[0] ) * bnorm;
				wnorm               += weights[m_Channels];
				wnorm                = 1.0f / wnorm;
				for( unsigned int c=0; c<m_Channels+1; ++c )
					weights[c] *= wnorm;
				// Now we have the weighting factors, we may compute X'WX; to reduce
				// the number of products, we compute first WX in an auxiliar matrix:
				LSMatrixType  WX( m_Channels+1 );
				for( unsigned int c=0; c<m_Channels+1; ++c ){
					for( unsigned int j=0; j<DTI_ESTIMATE_DIM; ++j )
						WX[c][j] = weights[c] * m_X[c][j];
				}
				// Now, we compute X'WX
				InverseLSType XWX;
				for( unsigned int i=0; i<DTI_ESTIMATE_DIM; ++i ){
					for( unsigned int j=i; j<DTI_ESTIMATE_DIM; ++j ){
						XWX[i][j] = itk::NumericTraits<double>::Zero;
						for( unsigned int c=0; c<m_Channels+1; ++c )
							XWX[i][j] += m_X[c][i] * WX[c][j];
						XWX[j][i] = XWX[i][j];
					}
				}
				// And the inverse of this matrix, (X'WX)^(-1)
				InverseLSType iXWX;
				if(   !this->ComputeInverseMatrix( XWX, iXWX )   ) // Matrix is singular; use the previous value
					break;
				// Once the inverse has been computed, we may compute the following
				// estimate of the components of the tensor; first, we compute
				// X'WY:
				for( unsigned int c=0; c<m_Channels+1; ++c ) // WY
					weights[c] *= dwi[c];
				double intermediate[DTI_ESTIMATE_DIM];     // X'WY
				for( unsigned int k=0; k<DTI_ESTIMATE_DIM; ++k ){
					intermediate[k] = itk::NumericTraits<double>::Zero;
					for( unsigned int c=0; c<m_Channels+1; ++c )
						intermediate[k] += m_X[c][k] * weights[c];
				}
				// At this point, we no longer need the vector of weights, so we
				// delete it:
				delete[] weights;
				// Finally, we compute: (X'WX)^(-1)(X'WY), which is, the solution to
				// the WLS problem:
				for( unsigned int k=0; k<DTI_ESTIMATE_DIM; ++k ){
					top[k] = itk::NumericTraits<double>::Zero;
					for( unsigned int c=0; c<DTI_ESTIMATE_DIM; ++c )
						top[k] += iXWX[k][c] * intermediate[c];
				}
				// The components of the tensor are now stored on op, and ready for
				// the next iteration
			}
			// Now, we have to undo the normalization to yield the final estimated tensor:
			double unorm = 10.0f/m_B;
			for( unsigned int k=0; k<DTI_ESTIMATE_DIM-1; ++k )
				op[k] = top[k+1]*unorm;
			// If we have to compute the T2 value, do it now!
			if( m_ComputeT2 ){
				t2.Set(   (float)::exp( 0.05f * top[0] )   );
				++t2;
			}

			// After a given number of iterations (maybe 0), we have an estimate of
			// the components of the tensor:
			it.Set( op );
		}
		else{
			op.Fill( itk::NumericTraits<double>::Zero );
			it.Set( op );
			if( m_ComputeT2 ){
				t2.Set( 0.0f );
				++t2;
			}
		}
		if( m_Mask != (MaskImagePointer)NULL ){
			++mit;
		}
	}
	// Delete auxiliar vector!!!
	delete[] dwi;
	return;
}
	
/**
 * The matrix is always 7x7 in size and regular, so we use Gaussian
 * elimination to invert it; any iterative method would require more
 * computation, since the number of iterations is large compared to the
 * size of the matrix. Besides, we may exploit that it is symetric.
 */
	
template< class TInputImage, class TOutputImage >
bool DTIEstimateTensorFilter< TInputImage, TOutputImage >
::ComputeInverseMatrix( InverseLSType input, InverseLSType& output ) const
{
	const double tol = 1e-5;
	// First of all, we initiallise the output matrix to be an identity matrix:
	output.Fill(  itk::NumericTraits<double>::Zero  );
	for( unsigned int col=0; col<DTI_ESTIMATE_DIM; ++col )
		output[col][col] = itk::NumericTraits<double>::One;
	// For each column col = 1 to 6, we need to make zeros in rows from
	// col+1 to 7 (note that in C++ array indices are 0-based):
	for( unsigned int col=0; col<DTI_ESTIMATE_DIM-1; ++col ){ // For each column
		// We need a non-null element in the position (col,col), in order to
		// accomplish gaussian elimination:
		if( fabs(input[col][col]) <= tol ){
			// Bad luck! The element is zero. We need to add a complete row to
			// the row in position c, so that the new element in position (c,c)
			// is not null. Find the first row for which the element (row,col) 
			// is non-zero:
			unsigned int row = col+1;
			while( fabs(input[row][col]) <= tol && row<DTI_ESTIMATE_DIM ){++row;}
			// If we are not able to find a row satisfying this condition, then
			// the matrix is singular, and this should not be the case; for
			// this reason, we do not perform bound checking, for efficiency. We
			// assume that row is a valid position, and then correct the input
			// and output:
			if( row==DTI_ESTIMATE_DIM ) // Singular matrix!!!
				return false;
			for( unsigned int cc=col; cc<DTI_ESTIMATE_DIM; ++cc )
				input[col][cc]  += input[row][cc];
			for( unsigned int cc=0; cc<DTI_ESTIMATE_DIM; ++cc )
				output[col][cc] += output[row][cc];
		}
		// At this point, we have a valid (col,col), element. We scale the whole
		// corresponding col-th row so that the pivoting element is simply 1:
		double scale = itk::NumericTraits<double>::One / input[col][col];
		for( unsigned int cc=col; cc<DTI_ESTIMATE_DIM; ++cc )
			input[col][cc]  *= scale;
		for( unsigned int cc=0; cc<DTI_ESTIMATE_DIM; ++cc )
			output[col][cc] *= scale;
		// Now, we may perform gaussian elimination for each row:
		for( unsigned int row=col+1; row<DTI_ESTIMATE_DIM; ++row ){ // For each row
			double scale = input[row][col]; // This is the scale, since input[col][col] = 1.
			// Once again, for each column, we add the corresponding scaled
			// version of the pivoting element; however, in the input matrix,
			// values at the left of this column are assumed to be already zero:
			for( unsigned int cc=col; cc<DTI_ESTIMATE_DIM; ++cc ) // Only the columns from col
				input[row][cc] -= scale * input[col][cc];
			for( unsigned int cc=0; cc<DTI_ESTIMATE_DIM; ++cc ) // All columns!
				output[row][cc] -= scale * output[col][cc];
			// We have completed this row
		}
		// We have completed this column
	}
	// Now we have an upper-triangular matrix, where all diagonal elements are
	// just 1, except for the last one; we need to normalise the last row:
	double scale = itk::NumericTraits<double>::One / input[DTI_ESTIMATE_DIM-1][DTI_ESTIMATE_DIM-1];
	input[DTI_ESTIMATE_DIM-1][DTI_ESTIMATE_DIM-1] = itk::NumericTraits<double>::One;
	for( unsigned int cc=0; cc<DTI_ESTIMATE_DIM; ++cc ) // All columns!
		output[DTI_ESTIMATE_DIM-1][cc] *= scale;
	// A similar process is now done to eliminate all elements above the
	// principal diagonal:
	// For each column col = 7 to 2, we need to make zeros in rows from
	// col-1 to 1 (note that in C++ array indices are 0-based):
	for( unsigned int col=DTI_ESTIMATE_DIM-1; col>0; --col ){ // For each column, reverse order
		// In this case, we KNOW that we already have a valid pivoting element,
		// which is simply 1, and we may perform gaussian elimination for 
		// each row:
		for( int row=col-1; row>=0; --row ){ // For each row, reverse order (row MUST be signed!!!)
			double scale = input[row][col];  // This is the scale, since input[col][col] = 1.
			// In this case, the pivoting row of the input matrix consists of
			// one single non-null element, which in fact is just "1", so in
			// each row of input we simply set the corresponding element to 0;
			// Indeed, since this value will not be used, we can even avoid
			// this assignment.
			// For the output matrix, we exploit the fact that it must be
			// symmetric, and compute only the elements below the main
			// diagonal.
			// This instruction can be droped, since input[row][col] is no longer
			// needed:
			// input[row][col] = itk::NumericTraits<double>::Zero;
			//for( unsigned int cc=0; cc<col; ++cc ){
			for( unsigned int cc=0; cc<=(unsigned int)row; ++cc )
				output[row][cc] -= scale * output[col][cc];
			// We have completed this row
		}
		// We have completed this column
	}
	// We exploit the symmetry to compute the values over the main diagonal,
	// since the inverted matrix MUST be symmetric:
	for( unsigned int i=0; i<DTI_ESTIMATE_DIM-1; ++i ){
		for( unsigned int j=i+1; j<DTI_ESTIMATE_DIM; ++j )
			output[i][j] = output[j][i];
	}
	// The matrix has been inverted!!!
	return true;
}

template< class TInputImage, class TOutputImage >
void DTIEstimateTensorFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData( void )
{

	// Precompute the inverse matrix for the LS problem
	if( m_Channels < DTI_ESTIMATE_DIM-1 )
		itkExceptionMacro( << "Not enough gradient directions to perform LS-fitting!!!" );
	// If this is the first call to Update(), the last equation, i.e., the equation corresponding
	// to the baseline image fit, has not been set. Check:
	if( m_X.size() <= m_Channels ){
		EquationType eq( itk::NumericTraits<double>::Zero );
		eq[0] = 0.05f;
		m_X.push_back( eq );
	}
	// Compute X'X
	InverseLSType xx;
	for( unsigned int i=0; i<DTI_ESTIMATE_DIM; ++i ){
		for( unsigned int j=i; j<DTI_ESTIMATE_DIM; ++j ){
			xx[i][j] = itk::NumericTraits<double>::Zero;
			for( unsigned int c=0; c<m_Channels+1; ++c )
				xx[i][j] += m_X[c][i] * m_X[c][j];
			xx[j][i] = xx[i][j];
		}
	}
	// Compute the inverse of X'X, (X'X)^(-1); it may sound absurd to pass an
	// attribute of the class as an argument to a method of this same class, but
	// this allows as to use this same method when it comes to compute the
	// inverse of the WLS problem at each voxel.
	if( !this->ComputeInverseMatrix( xx, m_Inverse ) )
		itkExceptionMacro( << "Matrix is too close to singular!!!" );
	// Finally, compute the pseudo inverse matrix of the problem, (X'X)^(-1)X'; 
	// the values of the components of the tensor are directly computed 
	// multiplying this matrix by the vecotr of unknowns. Note that this matrix 
	// is stored transposed.
	m_PseudoInverse = LSMatrixType( m_Channels+1 );
	for( unsigned int i=0; i<DTI_ESTIMATE_DIM; ++i ){
		for( unsigned int j=0; j<m_Channels+1; ++j ){
			m_PseudoInverse[j][i] = itk::NumericTraits<double>::Zero;
			for( unsigned int s=0; s<DTI_ESTIMATE_DIM; ++s )
				m_PseudoInverse[j][i] += m_Inverse[i][s] * m_X[j][s];
		}
	}
	if( m_ComputeT2 ){
		m_T2 = T2ImageType::New();
		m_T2->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
		m_T2->SetOrigin( this->GetInput()->GetOrigin() );
		m_T2->SetSpacing( this->GetInput()->GetSpacing() );
		m_T2->SetDirection( this->GetInput()->GetDirection() );
		m_T2->Allocate();
		m_T2->FillBuffer( 0.0f );
	}
	return;
}
	
	
} // end namespace itk


#endif
