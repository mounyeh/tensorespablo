/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineTransform.txx,v $
  Language:  C++
  Date:      $Date: 2005/10/05 18:20 $
  Version:   $Revision: 1.0 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

  For efficiency, no copy of parameters set by the setParameters() method is
  stored. I only keep a pointer to it and I assume that it is externally 
  mantained. After calling RefineGrid() the number of parameters needed is
  increased, so user must ensure that the externally mantained parameters 
  grow up to the new size. For this purpose, he may use the following sentence:
      
	  parameters.SetSize(transform->GetNumberOfParameters());

=========================================================================*/
#ifndef _itkBSplineTransform_txx
#define _itkBSplineTransform_txx

#include "itkBSplineTransform.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{

// Constructor with default arguments [10-5-2005-10:10]
template<class TScalarType, unsigned int NDimensions>
BSplineTransform<TScalarType, NDimensions>::BSplineTransform():Superclass(SpaceDimension,ParametersDimension)
{
  // This transform only suport 2D images:
  if(NDimensions != 2)
	  itkExceptionMacro( << "This transform is implemented only for 2 dimensions");
  // Point spline coefficients matrices to NULL value:
  PHIx = NULL;
  PHIy = NULL;
  m_ReturnedParameters = NULL;
  // Number of points in each direction:
  NX = 0;
  NY = 0;
  // Last visited grid point for computation of Jacobian:
  m_LastVisited[0] = 0;
  m_LastVisited[1] = 0;
}
    

// Destructor [10-5-2005-10:10]
template<class TScalarType, unsigned int NDimensions>
BSplineTransform<TScalarType, NDimensions>::~BSplineTransform()
{
    // Free allocated buffers:
	if((PHIx != NULL) && (PHIy != NULL)){
		for(unsigned int k=0; k<NX; k++){
			delete PHIx[k];
			delete PHIy[k];
		}
		delete(PHIx);
		delete(PHIy);
	}
	return;
}


// Set the number of points in the grid [10-5-2005-10:50]
template <class TScalarType, unsigned int NDimensions>
void BSplineTransform<TScalarType, NDimensions>::SetNumberOfPoints(const unsigned int pointsx, const unsigned int pointsy)
{
	// You must call this method before "SetParameters()"
	/* With N interior points of the grid (that is, we divide the width
	of the image by N), we need N+3 grid points. The points[0] value refers
	to the value N rather than the number of grid points, so it must be 
	interpreted as the factor we divide the width of the image by.*/
	// Update the value of the size of the grid:
	unsigned int k;
	
	// Free the old buffers, if needed:
	if((PHIx != NULL) && (PHIy != NULL)){
		for(k=0; k<NX; k++){
			delete PHIx[k];
			delete PHIy[k];
		}
		delete PHIx;
		delete PHIy;
	}
	
	NX = pointsx + 3;
	NY = pointsy + 3;
	
	// Allocate a buffer for the spline parameters:
	PHIx = new double*[NX];
	PHIy = new double*[NX];
	for(unsigned int k=0; k<NX; k++){
		PHIx[k] = new double[NY];
		PHIy[k] = new double[NY];
		// And initialize values to zero:
		for(unsigned int l=0; l<NY; l++){
			PHIx[k][l] = 0.0;
			PHIy[k][l] = 0.0;
		}
	}
	// Allocate buffer for Jacobian:
	m_Jacobian.set_size(2,2*NX*NY);
	m_Jacobian.Fill(NumericTraits<typename JacobianType::ValueType>::Zero);

	return;
}


// Set the region the grid extends over [10-5-2005-17:00]
template <class TScalarType, unsigned int NDimensions>
void BSplineTransform<TScalarType, NDimensions>::SetGridRegion(InputPointType origin, InputVectorType spacing, SizeType size)
{
	// You must call this method before "SetParameters()"
	/* Input parameters origin, spacing, and size can be obtained from the
	requested region to register via the corresponding "Get" methods */
	m_Origin[0] = (double)(origin[0]);
	m_Origin[1] = (double)(origin[1]);
	/* Note that the last point on each dimension does not match the
	corresponding m_Extent value, but it match the value m_Extent-spacing.
	It is purposely done in order to avoid out-of-bound problems. */
	m_Extent[0] = ((double)(size[0]))*((double)(spacing[0]));
	m_Extent[1] = ((double)(size[1]))*((double)(spacing[1]));
}


//Refine the grid to the next level of resolution [10-5-2005-11:50]
template <class TScalarType, unsigned int NDimensions>
void BSplineTransform<TScalarType, NDimensions>::RefineGrid(void)
{
	unsigned int p;
	// Allocate a new buffer for the spline parameters:
	double** PHI2x = new double*[2*NX-3];
	double** PHI2y = new double*[2*NX-3];
	for(p=0; p<2*NX-3; p++)
	{
		PHI2x[p] = new double[2*NY-3];
		PHI2y[p] = new double[2*NY-3];
	}

	/* Update the values of the new buffer in order to obtain an equivalent spline. 
	C++ 2d-array index (0,0) represents position (-1,-1) in the matrix of spline
	coefficients, so odd indices of 2D array represent even indices of the
	grid, and vice-versa */

	// Loop for even-even indices of the matrix of coefficients:
	for(unsigned int k=1; k<2*NX-3; k+=2){
		for(unsigned int l=1; l<2*NY-3; l+=2){
			// x-coordinate:
			PHI2x[k][l]  = ( PHIx[k/2][l/2] + PHIx[k/2][l/2+2] + PHIx[k/2+2][l/2] + PHIx[k/2+2][l/2+2] )/64.0;
			PHI2x[k][l] += 6.0*( PHIx[k/2][l/2+1] + PHIx[k/2+1][l/2] + PHIx[k/2+1][l/2+2] + PHIx[k/2+2][l/2+1] )/64.0;
			PHI2x[k][l] += 36.0*PHIx[k/2+1][l/2+1]/64.0;
			// y-coordinate:
			PHI2y[k][l]  = ( PHIy[k/2][l/2] + PHIy[k/2][l/2+2] + PHIy[k/2+2][l/2] + PHIy[k/2+2][l/2+2] )/64.0;
			PHI2y[k][l] += 6.0*( PHIy[k/2][l/2+1] + PHIy[k/2+1][l/2] + PHIy[k/2+1][l/2+2] + PHIy[k/2+2][l/2+1] )/64.0;
			PHI2y[k][l] += 36.0*PHIy[k/2+1][l/2+1]/64.0;
		}
	}
	// Loop for even-odd indices of the matrix of coefficients:
	for(unsigned int k=1; k<2*NX-3; k+=2){
		for(unsigned int l=0; l<2*NY-3; l+=2){
			// x-coordinate:
			PHI2x[k][l]  = ( PHIx[k/2][l/2] + PHIx[k/2][l/2+1] + PHIx[k/2+2][l/2] + PHIx[k/2+2][l/2+1] )/16.0;
			PHI2x[k][l] += 6.0*( PHIx[k/2+1][l/2] + PHIx[k/2+1][l/2+1] )/16.0;
			// y-coordinate:
			PHI2y[k][l]  = ( PHIy[k/2][l/2] + PHIy[k/2][l/2+1] + PHIy[k/2+2][l/2] + PHIy[k/2+2][l/2+1] )/16.0;
			PHI2y[k][l] += 6.0*( PHIy[k/2+1][l/2] + PHIy[k/2+1][l/2+1] )/16.0;
		}
	}
	// Loop for odd-even indices of the matrix of coefficients:
	for(unsigned int k=0; k<2*NX-3; k+=2){
		for(unsigned int l=1; l<2*NY-3; l+=2){
			// x-coordinate:
			PHI2x[k][l]  = ( PHIx[k/2][l/2] + PHIx[k/2][l/2+2] + PHIx[k/2+1][l/2] + PHIx[k/2+1][l/2+2] )/16.0;
			PHI2x[k][l] += 6.0*( PHIx[k/2][l/2+1] + PHIx[k/2+1][l/2+1] )/16.0;
			// y-coordinate:
			PHI2y[k][l]  = ( PHIy[k/2][l/2] + PHIy[k/2][l/2+2] + PHIy[k/2+1][l/2] + PHIy[k/2+1][l/2+2] )/16.0;
			PHI2y[k][l] += 6.0*( PHIy[k/2][l/2+1] + PHIy[k/2+1][l/2+1] )/16.0;
		}
	}
	// Loop for odd-odd indices of the matrix of coefficients:
	for(unsigned int k=0; k<2*NX-3; k+=2){
		for(unsigned int l=0; l<2*NY-3; l+=2){
			// x-coordinate.
			PHI2x[k][l]  = ( PHIx[k/2][l/2] + PHIx[k/2][l/2+1] + PHIx[k/2+1][l/2] + PHIx[k/2+1][l/2+1] )/4.0; 
			// y-coordinate:
			PHI2y[k][l]  = ( PHIy[k/2][l/2] + PHIy[k/2][l/2+1] + PHIy[k/2+1][l/2] + PHIy[k/2+1][l/2+1] )/4.0;
		}
	}
	// Free memory asociated with old buffer
	for(p=0; p<NX; p++)
	{
		delete PHIx[p];
		delete PHIy[p];
	}
	delete PHIx;
	delete PHIy;
	// Update the buffer of spline coefficients:
	PHIx = PHI2x;
	PHIy = PHI2y;
	// Update the buffer for Jacobian:
	m_Jacobian.set_size(2,2*(2*NX-3)*(2*NY-3));
	m_Jacobian.Fill(NumericTraits<typename JacobianType::ValueType>::Zero);
	// Update the values of NX and NY to the next level of resolution:
	NX = 2*NX - 3;
	NY = 2*NY - 3;
	
	return;
}


/** Interpolate given points */
template <class TScalarType, unsigned int NDimensions>
void BSplineTransform<TScalarType, NDimensions>
::InterpolatePoints( InputVector points, OutputVector values)
{
	/* This method assumes you want to interpolate given displacements 'values' at given
	locations 'points', and for that purpose you want to correct the existing PHI matrices
	(PHIx and PHIy are not reset). So, you can implement a multirresolution interpolation 
	algorithm in the manner:

	transform->SetGridRegion(···);
	transform->SetNumberOfPoints(NX,NY);

	transform->InterpolatePoints(points,values);
	transform->RefineGrid();

	transform->InterpolatePoints(points,values);
	transform->RefineGrid();

	...
	
	transform->InterpolatePoints(points,values);
	transform->RefineGrid();

	and the resut will not be the same if you do MX=NX=1 and then refine the grid three times
	(the final values will be NX=NY=8) as in the case you do MX=NX=4 and refine the grid once
	(the final value will be also NX=NY=8). With the second procedure, the entire process will
	be slightly slower but it will result in a much more soft deformation, with a completelly
	equivalent final error.
	*/
	//-----------------------------------------------------------------------------------------------
	// Counters:
	unsigned int k;
	unsigned int l;
	// Number of points for interpolation:
	unsigned int NP = points.size();
	// The PHI coefficients:
	double PHI[4][4];
	//-----------------------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------------------
	// The divider coefficients:
	double** OMEGA  = new double*[NX];
	// The weighting factors:
	double** DELTAx = new double*[NX]; 
	double** DELTAy = new double*[NX];
	// Zero the dividers and weighting factors:
	for(k=0; k<NX; k++){
		OMEGA[k] = new double[NY]; DELTAx[k] = new double[NY]; DELTAy[k] = new double[NY];
		for(l=0; l<NY; l++){
			OMEGA[k][l] = 0.0; DELTAx[k][l] = 0.0; DELTAy[k][l] = 0.0;
		}
	}
	//-----------------------------------------------------------------------------------------------
	

	//-----------------------------------------------------------------------------------------------
	// Current Interpolating point and current Interpolating value:
	InputPointType  cpoint;
	OutputPointType cvalue;
	// Grid indices:
	unsigned int i;
	unsigned int j;
	// s and t parameters:
	double x;
	double y;
	// wkl's sum:
	double wsq;
	//-----------------------------------------------------------------------------------------------


	//-----------------------------------------------------------------------------------------------
	// For each point (x_c,y_c,z_c) in P do:
	for(unsigned int p=0; p<NP; p++)
	{
		// Get the current point:
		cpoint = points[p];
		cvalue = (values[p]-(TransformPoint(cpoint)-static_cast<OutputPointType>(cpoint)));

		// Obtain grid coordinates and s and t parameters:
		//      1.- Transform InputType to double and normalize to grid coordinates:
		x = ((double)(NX-3))*( ( (double)(cpoint[0]) - m_Origin[0] ) / m_Extent[0] );
		y = ((double)(NY-3))*( ( (double)(cpoint[1]) - m_Origin[1] ) / m_Extent[1] );
		
		// If the point is outside the region the grid is defined on, we ignore that point:
		if( !(x<0 || x>=NX-3 || y<0 || y>=NY-3) )
		{
			//  2.- Get the i and j indices:
			i = (unsigned int)(floor(x));
			j = (unsigned int)(floor(y));
			//  3.- Get the s and t parameters:
			x = x - floor(x);
			y = y - floor(y);

			// Compute \phi_{kl}: and add it to the corresponding \delta and \omega locations
			PHI[0][0] = B0(x)*B0(y); PHI[0][1] = B0(x)*B1(y); PHI[0][2] = B0(x)*B2(y); PHI[0][3] = B0(x)*B3(y);
			PHI[1][0] = B1(x)*B0(y); PHI[1][1] = B1(x)*B1(y); PHI[1][2] = B1(x)*B2(y); PHI[1][3] = B1(x)*B3(y);
			PHI[2][0] = B2(x)*B0(y); PHI[2][1] = B2(x)*B1(y); PHI[2][2] = B2(x)*B2(y); PHI[2][3] = B2(x)*B3(y);
			PHI[3][0] = B3(x)*B0(y); PHI[3][1] = B3(x)*B1(y); PHI[3][2] = B3(x)*B2(y); PHI[3][3] = B3(x)*B3(y);
			// Compute wkl's sum:
			wsq = 0.0;
			for(k=0; k<4; k++){
				for(l=0; l<4; l++)
					wsq += (PHI[k][l])*(PHI[k][l]);
			}
			// Compute \omega and \delta:
			for(k=0; k<4; k++)
			{
				for(l=0; l<4; l++)
				{
					DELTAx[i+k][j+l] += (PHI[k][l])*(PHI[k][l])*(((PHI[k][l])*((double)(cvalue[0])))/wsq);
					DELTAy[i+k][j+l] += (PHI[k][l])*(PHI[k][l])*(((PHI[k][l])*((double)(cvalue[1])))/wsq);
					OMEGA[i+k][j+l]  += (PHI[k][l])*(PHI[k][l]);
				}
			}
		}
	}
	//-----------------------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------------------
	// For all i,j do:
	for(k=0; k<NX; k++){
		for(l=0; l<NY; l++){
			if( OMEGA[k][l] > 1e-12 ){
				PHIx[k][l] += (DELTAx[k][l]) / (OMEGA[k][l]);
				PHIy[k][l] += (DELTAy[k][l]) / (OMEGA[k][l]);
			}
		}
	}
	//-----------------------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------------------
	for(k=0; k<NX; k++){
		delete[] OMEGA[k]; delete[] DELTAx[k]; delete DELTAy[k];
	}
	delete[] OMEGA; delete[] DELTAx; delete DELTAy;
	//-----------------------------------------------------------------------------------------------

	return;
}


/** Interpolate given displacement image */
template <class TScalarType, unsigned int NDimensions>
void BSplineTransform<TScalarType, NDimensions>::InterpolateImage( InterpolatedImageType* image )
{
	//-----------------------------------------------------------------------------------------------
	// Counters:
	unsigned int k;
	unsigned int l;
	// The PHI coefficients:
	double PHI[4][4];
	//-----------------------------------------------------------------------------------------------


	//-----------------------------------------------------------------------------------------------
	// The divider coefficients:
	double** OMEGA  = new double*[NX];
	// The weighting factors:
	double** DELTAx = new double*[NX]; 
	double** DELTAy = new double*[NX];
	// Zero the dividers and weighting factors:
	for(k=0; k<NX; k++){
		OMEGA[k] = new double[NY]; DELTAx[k] = new double[NY]; DELTAy[k] = new double[NY];
		for(l=0; l<NY; l++){
			OMEGA[k][l] = 0.0; DELTAx[k][l] = 0.0; DELTAy[k][l] = 0.0;
		}
	}
	//-----------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------
	// Current Interpolating point and current Interpolating value:
	typename InterpolatedImageType::PointType cpoint;
	typename InterpolatedImageType::PixelType cvalue;
	OutputPointType op;
	// Grid indices:
	unsigned int i;
	unsigned int j;
	// s and t parameters:
	double x;
	double y;
	// wkl's sum:
	double wsq;
	//-----------------------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------------------
	// Image iterator:
	itk::ImageRegionIteratorWithIndex< InterpolatedImageType > it( image, image->GetRequestedRegion() );
	it.GoToBegin();

	// Covert from index to phisical point:
	typename InterpolatedImageType::PointType origin      = image->GetOrigin();
	SpacingType spacing                                   = image->GetSpacing();
	//-----------------------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------------------
	// For each point (x_c,y_c,z_c) do:
	while( ! it.IsAtEnd() )
	{
		// Get the current point:
		typename InterpolatedImageType::IndexType cindex = it.GetIndex();
		for( unsigned int p=0; p<InterpolatedImageType::ImageDimension; p++ )
			cpoint[p] = origin[p] + (spacing[p])*((double)(cindex[p]));
		cvalue = it.Get();
						
		op     = this->TransformPoint( cpoint );
		
		// Interpolate the error commited in the interpolation:
		for( unsigned int p=0; p<InterpolatedImageType::ImageDimension; p++ ){
			cvalue[p] += ( op[p] - cpoint[p] );
		}
			
		// Obtain grid coordinates and s and t parameters:
		//      1.- Transform InputType to double and normalize to grid coordinates:
		x = ((double)(NX-3))*( ( (double)(cpoint[0]) - m_Origin[0] ) / m_Extent[0] );
		y = ((double)(NY-3))*( ( (double)(cpoint[1]) - m_Origin[1] ) / m_Extent[1] );
		
		// If the point is outside the region the grid is defined on, we ignore that point:
		if( !(x<0 || x>=NX-3 || y<0 || y>=NY-3) )
		{
			//  2.- Get the i and j indices:
			i = (unsigned int)(floor(x));
			j = (unsigned int)(floor(y));
			//  3.- Get the s and t parameters:
			x = x - floor(x);
			y = y - floor(y);

			// Compute \phi_{kl}: and add it to the corresponding \delta and \omega locations
			PHI[0][0] = B0(x)*B0(y); PHI[0][1] = B0(x)*B1(y); PHI[0][2] = B0(x)*B2(y); PHI[0][3] = B0(x)*B3(y);
			PHI[1][0] = B1(x)*B0(y); PHI[1][1] = B1(x)*B1(y); PHI[1][2] = B1(x)*B2(y); PHI[1][3] = B1(x)*B3(y);
			PHI[2][0] = B2(x)*B0(y); PHI[2][1] = B2(x)*B1(y); PHI[2][2] = B2(x)*B2(y); PHI[2][3] = B2(x)*B3(y);
			PHI[3][0] = B3(x)*B0(y); PHI[3][1] = B3(x)*B1(y); PHI[3][2] = B3(x)*B2(y); PHI[3][3] = B3(x)*B3(y);
			// Compute wkl's sum:
			wsq = 0.0;
			for(k=0; k<4; k++){
				for(l=0; l<4; l++)
					wsq += (PHI[k][l])*(PHI[k][l]);
			}
			// Compute \omega and \delta:
			for(k=0; k<4; k++)
			{
				for(l=0; l<4; l++)
				{
					DELTAx[i+k][j+l] += (PHI[k][l])*(PHI[k][l])*(((PHI[k][l])*((double)(cvalue[0])))/wsq);
					DELTAy[i+k][j+l] += (PHI[k][l])*(PHI[k][l])*(((PHI[k][l])*((double)(cvalue[1])))/wsq);
					OMEGA[i+k][j+l]  += (PHI[k][l])*(PHI[k][l]);
				}
			}
		}
		++it;
	}
	//-----------------------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------------------
	// For all i,j do:
	for(k=0; k<NX; k++){
		for(l=0; l<NY; l++){
			if( OMEGA[k][l] > 1e-12 ){
				PHIx[k][l] -= (DELTAx[k][l]) / (OMEGA[k][l]);
				PHIy[k][l] -= (DELTAy[k][l]) / (OMEGA[k][l]);
			}
		}
	}
	//-----------------------------------------------------------------------------------------------


	//-----------------------------------------------------------------------------------------------
	for(k=0; k<NX; k++){
		delete[] OMEGA[k]; delete[] DELTAx[k]; delete DELTAy[k];
	}
	delete[] OMEGA; delete[] DELTAx; delete DELTAy;
	//-----------------------------------------------------------------------------------------------

	return;
}


// Get the parameters [10-5-2005-12:00]
template <class TScalarType, unsigned int NDimensions>
const typename BSplineTransform<TScalarType, NDimensions>::ParametersType &
BSplineTransform<TScalarType, NDimensions>::GetParameters( void ) const
{
	// Appropriatelly resize parameters buffer:
	(*m_ReturnedParameters).SetSize(2*NX*NY);
	// x-coordinate spline coefficients:
	for( unsigned int k=0; k<NX; k++ ){
		for( unsigned int l=0; l<NY; l++ ){
			(*m_ReturnedParameters)[NX*l + k] = PHIx[k][l];
		}
	}
	// y-coordinate spline coefficients:
	for( unsigned int k=0; k<NX; k++ ){
		for( unsigned int l=0; l<NY; l++ ){
			(*m_ReturnedParameters)[NX*NY + NX*l + k] = PHIy[k][l];
		}
	}
	return (*m_ReturnedParameters);
}


// Set the parameters [10-5-2005-12:30]
template <class TScalarType, unsigned int NDimensions>
void
BSplineTransform<TScalarType, NDimensions>::SetParameters( const ParametersType & parameters )
{
	/* As in BSplineDeformableTransform, we keep only a pointer to
	the externally managed parameters, but not the actual parameters, 
	because they are already in the PHI matrices*/
	
	if(parameters.Size() != 2*NX*NY){
		itkExceptionMacro(<< "Unexpected number of parameters. Needed: " << 2*NX*NY << 
		               ", provided: " << parameters.Size() << ". Make sure you have called SetNumberOfPoints().");
	}
	
	m_ReturnedParameters = (itk::Array<double>*)(&parameters);	
	
	for( unsigned int k=0; k<NX; k++ ){
		for( unsigned int l=0; l<NY; l++ ){
			PHIx[k][l] = parameters[NX*l + k];         // x-coordinate
			PHIy[k][l] = parameters[NX*NY + NX*l + k]; // y-coordinate
		}
	}	
	
	return;
}


// Set the parameters so the spline represent an affine transform [18-5-2005-17:15]
template <class TScalarType, unsigned int NDimensions>
void
BSplineTransform<TScalarType, NDimensions>::SetAffineTransform( const ParametersType & parameters )
{
	// In this method, we have to take into account that with 1x1 grid we normalize the input coordinates
	// so they fall into the interval [0,1), what leads to an extra transformation which makes us to modify
	// the base parameters
	
	if( (!NX) || (!NY))
		itkExceptionMacro( << "Please, initialize transform" );
	if(parameters.Size() != 6)
		itkExceptionMacro( << "You must provide 6 parameters for an affine transform");
	
	// Get the affine transform parameters:
	double l11 = parameters[0];
	double l12 = parameters[1];
	double l21 = parameters[2];
	double l22 = parameters[3];
	double Tx  = parameters[4];
	double Ty  = parameters[5];
	
	double PHIxBase = (l11-1.0)*m_Origin[0] + l12*m_Origin[1] + Tx;
	double PHIyBase = l21*m_Origin[0] + (l22-1.0)*m_Origin[1] + Ty;

	for (unsigned int k=0; k<NX; k++){
		for (unsigned int l=0; l<NY; l++){
			PHIx[k][l] = PHIxBase + ((double)k-1.0)*( m_Extent[0]/((double)NX-3.0)*(l11-1) )
				+ ((double)l-1.0)*( m_Extent[1]/((double)NY-3.0) )*l12;
			PHIy[k][l] = PHIyBase + ((double)l-1.0)*( m_Extent[1]/((double)NY-3.0)*(l22-1) )
				+ ((double)k-1.0)*( m_Extent[0]/((double)NX-3.0) )*l21;
		}
	}
	
	return;
}


// Get the parameters [10-5-2005-12:00]
template <class TScalarType, unsigned int NDimensions>
const typename BSplineTransform<TScalarType, NDimensions>::ParametersType
BSplineTransform<TScalarType, NDimensions>::GetInverseParameters( void ) const
{
	//----------------------------------------------------------------------------
	// Current parameters vector:
	ParametersType currentParameters = this->GetParameters();
	// Inverse parameters matrices:
	double** PSIx = new double*[NX];
	double** PSIy = new double*[NX];
	// Initialize parameters:
	for( unsigned int k=0; k<NX; k++ ){
		PSIx[k] = new double[NY]; PSIy[k] = new double[NY];
		for( unsigned int l=0; l<NY; l++ ){
			PSIx[k][l] = -PHIx[k][l]; PSIy[k][l] = -PHIy[k][l];
		}
	}
	//----------------------------------------------------------------------------


	//----------------------------------------------------------------------------
	// Original point and transformed point:
	InputPointType  point;
	InputPointType  point1;
	OutputPointType point2;
	// Iteratively find the inverse transformation at the grid points; we only consider
	// inner points, and assume that outside the deformation is constant in x and y.
	// Parameters for iterations:
	double tol        = 1e-5;
	double dif        = 2.0*tol;
	double cdif       = 0.0;
	unsigned int cont = 0;
	unsigned int MAX  = 10;
	//----------------------------------------------------------------------------


	//----------------------------------------------------------------------------
	// Coordinates of the grid point:
	while( dif>tol ){
		//----------------------------------------------------------------------------
		if( cont==MAX )
			break;
		else
			cont++;
		//----------------------------------------------------------------------------
		dif = 0.0;
		for( unsigned int k=1; k<NX-2; k++ ){
			for( unsigned int l=1; l<NY-2; l++ ){
				// Get the real-world coordinates of the grid point and transform them with the current iteration:
				point[0] = ((double)(k-1))/((double)(NX-3))*(m_Extent[0]) + m_Origin[0];
				point[1] = ((double)(l-1))/((double)(NY-3))*(m_Extent[1]) + m_Origin[1];
				point1[0] = point[0] + PSIx[k][l];
				point1[1] = point[1] + PSIy[k][l];
				// Neumann boundary conditions:
				if( point1[0]<m_Origin[0] ){point1[0] = m_Origin[0];}
				if( point1[0]>m_Origin[0]+m_Extent[0] ){point1[0] = m_Origin[0]+m_Extent[0];}
				if( point1[1]<m_Origin[1] ){point1[1] = m_Origin[1];}
				if( point1[1]>m_Origin[1]+m_Extent[1] ){point1[1] = m_Origin[1]+m_Extent[1];}
				// Transformation of the point with the original B-spline deformation:
				point2 = this->TransformPoint( point1 );
				// Current error:
				cdif =  (point[0]-point2[0])*(point[0]-point2[0])  +  (point[1]-point2[1])*(point[1]-point2[1]);
				if( cdif>dif )
					dif = cdif;
				// Fixed point updating rule:
				PSIx[k][l] = point1[0]-point2[0];
				PSIy[k][l] = point1[1]-point2[1];
			}
		}
	}
	//----------------------------------------------------------------------------


	//----------------------------------------------------------------------------
	// Extension to the borders:
	for( unsigned int k=0; k<NX; k++ ){
		PSIx[k][0] = PSIx[k][1]; PSIx[k][NY-1] = PSIx[k][NY-2];
		PSIy[k][0] = PSIy[k][1]; PSIy[k][NY-1] = PSIy[k][NY-2];
	}
	for( unsigned int k=0; k<NY; k++ ){
		PSIx[0][k] = PSIx[1][k]; PSIx[NX-1][k] = PSIx[NX-2][k];
		PSIy[0][k] = PSIy[1][k]; PSIy[NX-1][k] = PSIy[NX-2][k];
	}
	//----------------------------------------------------------------------------


	//----------------------------------------------------------------------------
	// Convert matrices into parameters:
	double* inverseXValues     = new double[NX*NY];
	double* inverseYValues     = new double[NX*NY];
	double* inverseXParameters = new double[NX*NY];
	double* inverseYParameters = new double[NX*NY];
	for( unsigned int k=0; k<NX; k++ ){
		for( unsigned int l=0; l<NY; l++ ){
			inverseXValues[NX*l + k]         = PSIx[k][l];
			inverseYValues[NX*l + k]         = PSIy[k][l];
		}
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	for( unsigned int k=0; k<NX; k++ ){
		delete[] PSIx[k]; delete[] PSIy[k];
	}
	delete[] PSIx; delete PSIy;
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// Solve the linear system to obtain grid weights from grid values:
	unsigned int sizes[2] = {NY,NX};
	this->Thomas( (double*)inverseXValues, (double*)inverseXParameters, (unsigned int*)sizes, 1 );
	this->Thomas( (double*)inverseYValues, (double*)inverseYParameters, (unsigned int*)sizes, 1 );
	delete[] inverseXValues;
	delete[] inverseYValues;
	//----------------------------------------------------------------------------


	//----------------------------------------------------------------------------
	// Return inverse parameters:
	ParametersType inverseParameters;
	inverseParameters.SetSize( this->GetNumberOfParameters() );
	for( unsigned int k=0; k<NX*NY; k++ ){
		inverseParameters[k]       = inverseXParameters[k];
		inverseParameters[k+NX*NY] = inverseYParameters[k];
	}
	delete[] inverseXParameters;
	delete[] inverseYParameters;
	//----------------------------------------------------------------------------


	return inverseParameters;
}

//Solve a block-tridiagonal system by means of the Thomas algorithm
template<class TScalarType, unsigned int NDimensions>
void BSplineTransform<TScalarType, NDimensions>::Thomas( double* values, double* results, unsigned int* sizes, unsigned int order ) const
{
	//----------------------------------------------------------------------------
	unsigned int DIM = 1;
	for( unsigned int k=0; k<order; k++ )
			DIM*=sizes[k];
	//----------------------------------------------------------------------------


	//----------------------------------------------------------------------------
	double*  gamma = new double[sizes[0]+1];
	double** theta = new double*[sizes[0]+1];
	for( unsigned int k=0; k<=sizes[0]; k++ )
		theta[k] = new double[DIM];
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// Forward swept to compute the recursive multiplicative factors:
	gamma[0] = 1;
	for( unsigned int d=0; d<DIM; d++ )
		theta[0][d] = 0.0;
	for( unsigned int k=0; k<sizes[0]; k++ ){
		gamma[k+1] = -1.0/( gamma[k] + 4.0 );
		if( order==0 )
			theta[k+1][0] = (-gamma[k+1])*( 6.0*values[k] - theta[k][0] );
		else{
			this->Thomas( ((double*)values)+(k*DIM), (double*)(theta[k+1]), (unsigned int*)(sizes+1), order-1 );
			for( unsigned int d=0; d<DIM; d++ )
				theta[k+1][d] = (-gamma[k+1])*( 6.0*theta[k+1][d] - theta[k][d] );
		}
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// Backward swept
	for( unsigned int d=0; d<DIM; d++ )
		results[(sizes[0]-1)*DIM+d] = (  (theta[sizes[0]][d])/(1.0-gamma[sizes[0]])  )*gamma[sizes[0]] + theta[sizes[0]][d];
	for( int k=sizes[0]-2; k>=0; k-- ){
		for( unsigned int d=0; d<DIM; d++ )
			results[k*DIM+d] = (gamma[k+1])*(results[(k+1)*DIM+d]) + (theta[k+1][d]);
	}
	//----------------------------------------------------------------------------


	//----------------------------------------------------------------------------
	for( unsigned int k=0; k<=sizes[0]; k++ )
		delete[] theta[k];
	delete[] theta;
	delete[] gamma;
	//----------------------------------------------------------------------------
}


// Transform a point [10-5-2005-17:40]
template<class TScalarType, unsigned int NDimensions>
typename BSplineTransform<TScalarType, NDimensions>::OutputPointType
BSplineTransform<TScalarType, NDimensions>::TransformPoint(const InputPointType &point) const 
{
	// Offset in both directions:
	InputVectorType offset;
	// Displacement along each coordinate:
	double xoffset = 0.0;
	double yoffset = 0.0;
	// Grid indices:
	unsigned int i;
	unsigned int j;
	// Conversion to double type:
	double x = (double)(point[0]);
	double y = (double)(point[1]);	
	
	// If the point is outside the region the grid is defined on, we return the same point:
	if( x<m_Origin[0] || x>m_Origin[0]+m_Extent[0] || y<m_Origin[1] || y>m_Origin[1]+m_Extent[1]){
		return point;
	}

	// Normalization of coordinates to grid coordinates:
	x = ((double)(NX-3))*( ( x - m_Origin[0] ) / m_Extent[0] );
	y = ((double)(NY-3))*( ( y - m_Origin[1] ) / m_Extent[1] );
	
	// Get the i and j indices:
	i = (unsigned int)(floor(x));
	j = (unsigned int)(floor(y));

	// Get the s and t parameters:
	x = x - floor(x);
	y = y - floor(y);

	// Avoid border artifacts:
	if( i>NX-4 ){ i = NX-4; x = 1.0; }
	if( j>NY-4 ){ j = NY-4; y = 1.0; }
	
	// Get the displacements by means of the evaluation of splines:
	xoffset += PHIx[i][j]*B0(x)*B0(y);   xoffset += PHIx[i][j+1]*B0(x)*B1(y);   xoffset += PHIx[i][j+2]*B0(x)*B2(y);   xoffset += PHIx[i][j+3]*B0(x)*B3(y);
	xoffset += PHIx[i+1][j]*B1(x)*B0(y); xoffset += PHIx[i+1][j+1]*B1(x)*B1(y); xoffset += PHIx[i+1][j+2]*B1(x)*B2(y); xoffset += PHIx[i+1][j+3]*B1(x)*B3(y);
	xoffset += PHIx[i+2][j]*B2(x)*B0(y); xoffset += PHIx[i+2][j+1]*B2(x)*B1(y); xoffset += PHIx[i+2][j+2]*B2(x)*B2(y); xoffset += PHIx[i+2][j+3]*B2(x)*B3(y);
	xoffset += PHIx[i+3][j]*B3(x)*B0(y); xoffset += PHIx[i+3][j+1]*B3(x)*B1(y); xoffset += PHIx[i+3][j+2]*B3(x)*B2(y); xoffset += PHIx[i+3][j+3]*B3(x)*B3(y);

	yoffset += PHIy[i][j]*B0(x)*B0(y);   yoffset += PHIy[i][j+1]*B0(x)*B1(y);   yoffset += PHIy[i][j+2]*B0(x)*B2(y);   yoffset += PHIy[i][j+3]*B0(x)*B3(y);
	yoffset += PHIy[i+1][j]*B1(x)*B0(y); yoffset += PHIy[i+1][j+1]*B1(x)*B1(y); yoffset += PHIy[i+1][j+2]*B1(x)*B2(y); yoffset += PHIy[i+1][j+3]*B1(x)*B3(y);
	yoffset += PHIy[i+2][j]*B2(x)*B0(y); yoffset += PHIy[i+2][j+1]*B2(x)*B1(y); yoffset += PHIy[i+2][j+2]*B2(x)*B2(y); yoffset += PHIy[i+2][j+3]*B2(x)*B3(y);
	yoffset += PHIy[i+3][j]*B3(x)*B0(y); yoffset += PHIy[i+3][j+1]*B3(x)*B1(y); yoffset += PHIy[i+3][j+2]*B3(x)*B2(y); yoffset += PHIy[i+3][j+3]*B3(x)*B3(y);
	
	// Update the value of the point:
	offset[0] = xoffset;
	offset[1] = yoffset;
	
	// return:
	return (point + offset);
}
  

// Compute the Jacobian in one position [10-5-2005-18:20]
template<class TScalarType, unsigned int NDimensions>
const typename BSplineTransform<TScalarType, NDimensions>::JacobianType & 
BSplineTransform< TScalarType, NDimensions >::GetJacobian( const InputPointType &point) const
{
	unsigned int k;
	unsigned int l;

	unsigned int i;
	unsigned int j;
	// Conversion to double type:
	double x = (double)(point[0]);
	double y = (double)(point[1]);
	
	// First, we set to 0 the values of Jacobian, when needed:
	for(k=m_LastVisited[0]; k<m_LastVisited[0]+4; k++){
		for(l=m_LastVisited[1]; l<m_LastVisited[1]+4; l++){
			m_Jacobian(0,k+l*NX)       = 0.0;
			m_Jacobian(1,k+l*NX+NX*NY) = 0.0;
		}
	}

	// If the point is outside the region the grid is defined on, we return the same point:
	if( x<m_Origin[0] || x>m_Origin[0]+m_Extent[0] || y<m_Origin[1] || y>m_Origin[1]+m_Extent[1])
		return m_Jacobian;

	// Normalization of coordinates to grid coordinates:
	x = ((double)(NX-3))*( ( x - m_Origin[0] ) / m_Extent[0] );
	y = ((double)(NY-3))*( ( y - m_Origin[1] ) / m_Extent[1] );
	
	// Get the i and j indices:
	i = (unsigned int)(floor(x));
	j = (unsigned int)(floor(y));
	
	// Get the s and t parameters:
	x = x - floor(x);
	y = y - floor(y);

	// Update the LastVisited:
	m_LastVisited[0] = i;
	m_LastVisited[1] = j;
	
	// Fill the rest of Jacobian values (32 values):
	m_Jacobian(1,i+j*NX+NX*NY)         = B0(x)*B0(y);     m_Jacobian(0,i+j*NX)         = B0(x)*B0(y);
	m_Jacobian(1,(i+1)+j*NX+NX*NY)     = B1(x)*B0(y);     m_Jacobian(0,(i+1)+j*NX)     = B1(x)*B0(y);
	m_Jacobian(1,(i+2)+j*NX+NX*NY)     = B2(x)*B0(y);     m_Jacobian(0,(i+2)+j*NX)     = B2(x)*B0(y);
	m_Jacobian(1,(i+3)+j*NX+NX*NY)     = B3(x)*B0(y);     m_Jacobian(0,(i+3)+j*NX)     = B3(x)*B0(y);

	m_Jacobian(1,i+(j+1)*NX+NX*NY)     = B0(x)*B1(y);     m_Jacobian(0,i+(j+1)*NX)     = B0(x)*B1(y);
	m_Jacobian(1,(i+1)+(j+1)*NX+NX*NY) = B1(x)*B1(y);     m_Jacobian(0,(i+1)+(j+1)*NX) = B1(x)*B1(y);
	m_Jacobian(1,(i+2)+(j+1)*NX+NX*NY) = B2(x)*B1(y);     m_Jacobian(0,(i+2)+(j+1)*NX) = B2(x)*B1(y);
	m_Jacobian(1,(i+3)+(j+1)*NX+NX*NY) = B3(x)*B1(y);     m_Jacobian(0,(i+3)+(j+1)*NX) = B3(x)*B1(y);

	m_Jacobian(1,i+(j+2)*NX+NX*NY)     = B0(x)*B2(y);     m_Jacobian(0,i+(j+2)*NX)     = B0(x)*B2(y);
	m_Jacobian(1,(i+1)+(j+2)*NX+NX*NY) = B1(x)*B2(y);     m_Jacobian(0,(i+1)+(j+2)*NX) = B1(x)*B2(y);
	m_Jacobian(1,(i+2)+(j+2)*NX+NX*NY) = B2(x)*B2(y);     m_Jacobian(0,(i+2)+(j+2)*NX) = B2(x)*B2(y);
	m_Jacobian(1,(i+3)+(j+2)*NX+NX*NY) = B3(x)*B2(y);     m_Jacobian(0,(i+3)+(j+2)*NX) = B3(x)*B2(y);

	m_Jacobian(1,i+(j+3)*NX+NX*NY)     = B0(x)*B3(y);     m_Jacobian(0,i+(j+3)*NX)     = B0(x)*B3(y);
	m_Jacobian(1,(i+1)+(j+3)*NX+NX*NY) = B1(x)*B3(y);     m_Jacobian(0,(i+1)+(j+3)*NX) = B1(x)*B3(y);
	m_Jacobian(1,(i+2)+(j+3)*NX+NX*NY) = B2(x)*B3(y);     m_Jacobian(0,(i+2)+(j+3)*NX) = B2(x)*B3(y);
	m_Jacobian(1,(i+3)+(j+3)*NX+NX*NY) = B3(x)*B3(y);     m_Jacobian(0,(i+3)+(j+3)*NX) = B3(x)*B3(y);

	// return:
	return m_Jacobian;
}


// Compute the penalty term
template<class TScalarType, unsigned int NDimensions>
double BSplineTransform< TScalarType, NDimensions >::GetCorrection()
{
	// Initialization ckecking:
	if( !(NX*NY) )
		itkExceptionMacro( << "You must initialize the transform with SetNumberOfPoints() before you call GetCorrection()");

	//------------------------------------------------------------------------------------------------------------------------------
	// I take advantage on the fact that all double integrals are separable, so they can be expressed as a product of an
	// integral in x by an integral in y, whose values can indeed be precomputed

	// Pre-computed integrals:
	double D0D0[4][4];    // Zero-order derivatives
	double D1D1[4][4];    // First derivatives
	double D2D2[4][4];    // Second derivatives

	// int(B[m](s)·B[n](s), s=0..1)
	D0D0[0][0]=1.0/252.0;   D0D0[0][1]=43.0/1680.0;  D0D0[0][2]=1.0/84.0;     D0D0[0][3]=1.0/5040.0;
	D0D0[1][0]=43.0/1680.0; D0D0[1][1]=33.0/140.0;   D0D0[1][2]=311.0/1680.0; D0D0[1][3]=1.0/84.0;
	D0D0[2][0]=1.0/84.0;    D0D0[2][1]=311.0/1680.0; D0D0[2][2]=33.0/140.0;   D0D0[2][3]=43.0/1680.0;
	D0D0[3][0]=1.0/5040.0;  D0D0[3][1]=1.0/84.0;     D0D0[3][2]=43.0/1680.0;  D0D0[3][3]=1.0/252.0;

	// int(B[m]'(s)·B[n]'(s), s=0..1)
	D1D1[0][0]=1.0/20.0;    D1D1[0][1]=7.0/120.0;    D1D1[0][2]=-1.0/10.0;    D1D1[0][3]=-1.0/120.0;
	D1D1[1][0]=7.0/120.0;   D1D1[1][1]=17.0/60.0;    D1D1[1][2]=-29.0/120.0;  D1D1[1][3]=-1.0/10.0;
	D1D1[2][0]=-1.0/10.0;   D1D1[2][1]=-29.0/120.0;  D1D1[2][2]=17.0/60.0;    D1D1[2][3]=7.0/120.0;
	D1D1[3][0]=-1.0/120.0;  D1D1[3][1]=-1.0/10.0;    D1D1[3][2]=7.0/120.0;    D1D1[3][3]=1.0/20.0;

	// int(B[m]''(s)·B[n]''(s), s=0..1)
	D2D2[0][0]=1.0/3.0;     D2D2[0][1]=-0.5;         D2D2[0][2]=0.0;          D2D2[0][3]=1.0/6.0;
	D2D2[1][0]=-0.5;        D2D2[1][1]=1.0;          D2D2[1][2]=-0.5;         D2D2[1][3]=0.0;
	D2D2[2][0]=0.0;         D2D2[2][1]=-0.5;         D2D2[2][2]=1.0;          D2D2[2][3]=-0.5;
	D2D2[3][0]=1.0/6.0;     D2D2[3][1]=0.0;          D2D2[3][2]=-0.5;         D2D2[3][3]=1.0/3.0;

	//------------------------------------------------------------------------------------------------------------------------------

	// Normalization constants:
	double nx  = (double)NX-3.0;
	double ny  = (double)NY-3.0;
	double X   = (double)(m_Extent[0]);
	double Y   = (double)(m_Extent[1]);
	double la1 = (nx*nx*nx)/(ny*X*X*X*X);
	double la2 = (ny*ny*ny)/(nx*Y*Y*Y*Y);
	double la3 = 2.0*(nx*ny)/(X*X*Y*Y);
	
	// Reckoning:
	double val = 0.0;
	for( unsigned int i=0; i<NX-3; i++ ){ // 1
		for ( unsigned int j=0; j<NY-3; j++ ){ // 22
			// For each parameter (i, j), compute the integral in the cell:
			for( unsigned int k=0; k<4; k++ ){ // 333
				for( unsigned int l=0; l<4; l++ ){ // 4444
					// First sums
					for( unsigned int p=0; p<=3; p++ ){ // 55555
						for( unsigned int q=0; q<=3; q++ ){ // 666666
							// Second sums
							val += ( (PHIx[i+k][j+l])*(PHIx[i+p][j+q]) + (PHIy[i+k][j+l])*(PHIy[i+p][j+q]) )*
								( la1*(D2D2[k][p])*(D0D0[l][q]) + la2*(D0D0[k][p])*(D2D2[l][q]) + la3*(D1D1[k][p])*(D1D1[l][q]) );
						} // 666666
					} // 55555
				} // 4444
			} // 333
		} // 22
	} // 1

	// return:
	return val;
}


// Compute the Jacobian of the penalty term[10-5-2005-18:20]
template<class TScalarType, unsigned int NDimensions>
typename BSplineTransform<TScalarType, NDimensions>::CorrectionJacobianType 
BSplineTransform< TScalarType, NDimensions >::GetCorrectionJacobian()
{
	// Initialization ckecking:
	if( !(NX*NY) )
		itkExceptionMacro( << "You must initialize the transform with SetNumberOfPoints() before you call GetCorrectionJacobian()");

	//------------------------------------------------------------------------------------------------------------------------------
	// I take advantage on the fact that all double integrals are separable, so they can be expressed as a product of an
	// integral in x by an integral in y, whose values can indeed be precomputed

	// Pre-computed integrals:
	double D0D0[4][4];    // Zero-order derivatives
	double D1D1[4][4];    // First derivatives
	double D2D2[4][4];    // Second derivatives

	// int(B[m](s)·B[n](s), s=0..1)
	D0D0[0][0]=1.0/252.0;   D0D0[0][1]=43.0/1680.0;  D0D0[0][2]=1.0/84.0;     D0D0[0][3]=1.0/5040.0;
	D0D0[1][0]=43.0/1680.0; D0D0[1][1]=33.0/140.0;   D0D0[1][2]=311.0/1680.0; D0D0[1][3]=1.0/84.0;
	D0D0[2][0]=1.0/84.0;    D0D0[2][1]=311.0/1680.0; D0D0[2][2]=33.0/140.0;   D0D0[2][3]=43.0/1680.0;
	D0D0[3][0]=1.0/5040.0;  D0D0[3][1]=1.0/84.0;     D0D0[3][2]=43.0/1680.0;  D0D0[3][3]=1.0/252.0;

	// int(B[m]'(s)·B[n]'(s), s=0..1)
	D1D1[0][0]=1.0/20.0;    D1D1[0][1]=7.0/120.0;    D1D1[0][2]=-1.0/10.0;    D1D1[0][3]=-1.0/120.0;
	D1D1[1][0]=7.0/120.0;   D1D1[1][1]=17.0/60.0;    D1D1[1][2]=-29.0/120.0;  D1D1[1][3]=-1.0/10.0;
	D1D1[2][0]=-1.0/10.0;   D1D1[2][1]=-29.0/120.0;  D1D1[2][2]=17.0/60.0;    D1D1[2][3]=7.0/120.0;
	D1D1[3][0]=-1.0/120.0;  D1D1[3][1]=-1.0/10.0;    D1D1[3][2]=7.0/120.0;    D1D1[3][3]=1.0/20.0;

	// int(B[m]''(s)·B[n]''(s), s=0..1)
	D2D2[0][0]=1.0/3.0;     D2D2[0][1]=-0.5;         D2D2[0][2]=0.0;          D2D2[0][3]=1.0/6.0;
	D2D2[1][0]=-0.5;        D2D2[1][1]=1.0;          D2D2[1][2]=-0.5;         D2D2[1][3]=0.0;
	D2D2[2][0]=0.0;         D2D2[2][1]=-0.5;         D2D2[2][2]=1.0;          D2D2[2][3]=-0.5;
	D2D2[3][0]=1.0/6.0;     D2D2[3][1]=0.0;          D2D2[3][2]=-0.5;         D2D2[3][3]=1.0/3.0;

	//------------------------------------------------------------------------------------------------------------------------------

	// Normalization constants:
	double nx  = (double)NX-3.0;
	double ny  = (double)NY-3.0;
	double X   = (double)(m_Extent[0]);
	double Y   = (double)(m_Extent[1]);
	double la1 = (nx*nx*nx)/(ny*X*X*X*X);
	double la2 = (ny*ny*ny)/(nx*Y*Y*Y*Y);
	double la3 = 2.0*(nx*ny)/(X*X*Y*Y);

	// Declare the Jacobian and set the appropriate size:
	CorrectionJacobianType cj; 
	cj.set_size(2*NX*NY);

	// Cumulative values of each component of Jacobian
	double valx;
	double valy;
	// Define the lower limit of the first sums:
	unsigned int k0;
	unsigned int l0;
	// Define the upper limit of the first sums:
	unsigned int k3;
	unsigned int l3;

	// Reckoning of each Jacobian value:
	k0 = 0; k3 = 0;
	for( unsigned int r=0; r<NX; r++ ){ // 1
		if( r>=NX-3 ) k0++;      // Update the lower limit of first sumations, if needed
		l0 = 0; l3 = 0;
		for ( unsigned int s=0; s<NY; s++ ){ // 22
			if( s>=NY-3 ) l0++;  // Update the lower limit of first sumations, if needed
			valx = 0.0;
			valy = 0.0;
			for( unsigned int k=k0; k<=k3; k++ ){ // 333
				for( unsigned int l=l0; l<=l3; l++ ){ // 4444
					// sumation over each influence cell
					for( unsigned int p=0; p<=3; p++ ){ // 55555
						for( unsigned int q=0; q<=3; q++ ){ // 666666
							// Pre-computed integral in each influence cell
							valx += 2.0*(PHIx[p+r-k][q+s-l])*( la1*(D2D2[k][p])*(D0D0[l][q]) + la2*(D0D0[k][p])*(D2D2[l][q]) + la3*(D1D1[k][p])*(D1D1[l][q]) );
							valy += 2.0*(PHIy[p+r-k][q+s-l])*( la1*(D2D2[k][p])*(D0D0[l][q]) + la2*(D0D0[k][p])*(D2D2[l][q]) + la3*(D1D1[k][p])*(D1D1[l][q]) );
						} // 666666
					} // 55555
				} // 4444
			} // 333
			//-----------------------------------------------------------------------------------------------
			// Once here, we've got already the value of both Jacobian values, except for the normalization constant:
			cj[NX*s+r]       = valx;
			cj[NX*NY+NX*s+r] = valy;
			//-----------------------------------------------------------------------------------------------
			if( s<3 ) l3++;  // Update the upper limit of first sumations, if needed
		} // 22
		if( r<3 ) k3++;      // Update the upper limit of first sumations, if needed
	} // 1

	// return:
	return cj;
}

// Get the total number of parameters needed [10-5-2005-9:40]
template<class TScalarType, unsigned int NDimensions>
unsigned int BSplineTransform<TScalarType, NDimensions>::GetNumberOfParameters(void) const
{
	return 2*NX*NY;
}

// Set the parameters for an Identity transform of this class [10-5-2005-10:10]
template<class TScalarType, unsigned int NDimensions>
void
BSplineTransform<TScalarType, NDimensions>::
SetIdentity()
{
	if(NX*NY){
		for(unsigned int k=0; k<NX; k++){
			for(unsigned int l=0; l<NY; l++){
				PHIx[k][l] = 0.0;
				PHIy[k][l] = 0.0;
			}
		}
	}
	return;
}

// Print self [10-5-2005-10:10]
template<class TScalarType, unsigned int NDimensions>
void
BSplineTransform<TScalarType, NDimensions>::PrintSelf(std::ostream &os, Indent indent) const 
{
  //Superclass::PrintSelf(os,indent);
  
  os << indent << "%NX: " << NX << std::endl;
  os << indent << "%NY: " << NY << std::endl;
  os << indent << "%m_Origin = [" << m_Origin[0] << ", " << m_Origin[1] << "]" << std::endl;
  os << indent << "%m_Extent = [" << m_Extent[0] << ", " << m_Extent[1] << "]" << std::endl;
  os << indent << "PHIx = [" << std::endl;
  for(unsigned int k=0; k<NX; k++){
	  os << indent;
	  for(unsigned int l=0; l<NY; l++)
		  os << PHIx[k][l] << ", ";
	  os << std::endl;
  }
  os << indent << "];" << std::endl;
  os << indent << "PHIy = [" << std::endl;
  for(unsigned int k=0; k<NX; k++){
	  os << indent;
	  for(unsigned int l=0; l<NY; l++)
		  os << PHIy[k][l] << ", ";
	  os << std::endl;
  }
  os << indent << "];" << std::endl;
}
 
  
} // namespace

#endif
