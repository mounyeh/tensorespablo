/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkD3BSplineTransform.txx,v $
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
#ifndef _itkD3BSplineTransform_txx
#define _itkD3BSplineTransform_txx

#include "itkD3BSplineTransform.h"


namespace itk
{

// Constructor with default arguments [20-5-2005-9:25]
template<class TScalarType, unsigned int NDimensions>
D3BSplineTransform<TScalarType, NDimensions>::
D3BSplineTransform():Superclass(SpaceDimension,ParametersDimension)
{
  // This transform only suport D3 images::
  if(NDimensions != 3)
	  itkExceptionMacro( << "This transform is implemented only for 3 dimensions");
  // Point spline coefficients matrices to NULL value:
  PHIx = NULL;
  PHIy = NULL;
  PHIz = NULL;
  m_ReturnedParameters = NULL;
  // Number of points in each direction:
  NX = 0;
  NY = 0;
  NZ = 0;
  // Last visited grid point for computation of Jacobian:
  m_LastVisited[0] = 0;
  m_LastVisited[1] = 0;
  m_LastVisited[2] = 0;
}
    

// Destructor [20-5-2005-9:30]
template<class TScalarType, unsigned int NDimensions>
D3BSplineTransform<TScalarType, NDimensions>::
~D3BSplineTransform()
{
    // Free allocated buffers:
	if((PHIx != NULL) && (PHIy != NULL)){
		for(unsigned int k=0; k<NX; k++){
			for(unsigned int l=0; l<NY; l++){
				delete PHIx[k][l];
				delete PHIy[k][l];
				delete PHIz[k][l];
			}
			delete PHIx[k];
			delete PHIy[k];
			delete PHIz[k];
		}
		delete PHIx;
		delete PHIy;
		delete PHIz;
	}
	return;
}





// Set the number of points in the grid [20-5-2005-9:40]
template <class TScalarType, unsigned int NDimensions>
void D3BSplineTransform<TScalarType, NDimensions>
::SetNumberOfPoints(const unsigned int pointsx, const unsigned int pointsy, const unsigned int pointsz)
{
	// You must call this method before "SetParameters()"
	/* With N interior points of the grid (that is, we divide the width
	of the image by N), we need N+3 grid points. The pointsx value refers
	to the value N rather than the number of grid points, so it must be 
	interpreted as the factor we divide the width of the image by.*/
	// Update the value of the size of the grid:
	unsigned int k;
	unsigned int l;
	unsigned int m;
	
	// Free the old buffers, if needed:
	if((PHIx != NULL) && (PHIy != NULL)){
		for(k=0; k<NX; k++){
			for(l=0; l<NY; l++){
				delete PHIx[k][l];
				delete PHIy[k][l];
				delete PHIz[k][l];
			}
			delete PHIx[k];
			delete PHIy[k];
			delete PHIz[k];
		}
		delete PHIx;
		delete PHIy;
		delete PHIz;
	}

	NX = pointsx + 3;
	NY = pointsy + 3;
	NZ = pointsz + 3;
	
	// Allocate a buffer for the spline parameters:
	PHIx = new double**[NX];
	PHIy = new double**[NX];
	PHIz = new double**[NX];
	for(k=0; k<NX; k++){
		PHIx[k] = new double*[NY];
		PHIy[k] = new double*[NY];
		PHIz[k] = new double*[NY];
		for(l=0; l<NY; l++){
			PHIx[k][l] = new double[NZ];
			PHIy[k][l] = new double[NZ];
			PHIz[k][l] = new double[NZ];
			for(m=0; m<NZ; m++){
				// And initialize values to zero:
				PHIx[k][l][m] = 0.0;
				PHIy[k][l][m] = 0.0;
				PHIz[k][l][m] = 0.0;
			}
		}
	}
	// Allocate buffer for Jacobian:
	m_Jacobian.set_size(3,3*NX*NY*NZ);
	m_Jacobian.Fill(NumericTraits<typename JacobianType::ValueType>::Zero);

	return;
}




// Set the region the grid extends over [20-5-2005-9:40]
template <class TScalarType, unsigned int NDimensions>
void D3BSplineTransform<TScalarType, NDimensions>
::SetGridRegion(InputPointType origin, InputVectorType spacing, SizeType size)
{
	// You must call this method before "SetParameters()"
	/* Input parameters origin, spacing, and size can be obtained from the
	requested region to register via the corresponding "Get" methods */
	m_Origin[0] = (double)(origin[0]);
	m_Origin[1] = (double)(origin[1]);
	m_Origin[2] = (double)(origin[2]);
	/* Note that the last point on each dimension does not match the
	corresponding m_Extent value, but it match the value m_Extent-spacing.
	It is purposely done in order to avoid out-of-bound problems. */
	m_Extent[0] = ((double)(size[0]))*((double)(spacing[0]));
	m_Extent[1] = ((double)(size[1]))*((double)(spacing[1]));
	m_Extent[2] = ((double)(size[2]))*((double)(spacing[2]));
}




//Refine the grid to the next level of resolution [23-5-2005-12:10]
template <class TScalarType, unsigned int NDimensions>
void D3BSplineTransform<TScalarType, NDimensions>
::RefineGrid(void)
{
	//-------------------------------------------------------------------------------------
	// Interpolation filters:
	double v[3];
	double f[2][3];
	f[0][0] = 0.000; f[0][1] = 0.500; f[0][2] = 0.500;
	f[1][0] = 0.125; f[1][1] = 0.750; f[1][2] = 0.125;
	// Auxiliar index to avoid out-of-bounds troubles:
	int index[2];
	index[0] = 0; index[1] = 1;
	//-------------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------------
	unsigned int p;
	unsigned int q;
	// Allocate a new buffer for the spline parameters:
	double*** PHI2x = new double**[2*NX-3];
	double*** PHI2y = new double**[2*NX-3];
	double*** PHI2z = new double**[2*NX-3];
	for(p=0; p<2*NX-3; p++)
	{
		PHI2x[p] = new double*[2*NY-3];
		PHI2y[p] = new double*[2*NY-3];
		PHI2z[p] = new double*[2*NY-3];
		for(q=0; q<2*NY-3; q++){
			PHI2x[p][q] = new double[2*NZ-3];
			PHI2y[p][q] = new double[2*NZ-3];
			PHI2z[p][q] = new double[2*NZ-3];
		}
	}
	//-------------------------------------------------------------------------------------




	//-------------------------------------------------------------------------------------
	// Auxiliar buffer for succesive convolutions:
	double PHIC1x[NX][NY][NZ]; double PHIC2x[NX][NY][NZ];
	double PHIC1y[NX][NY][NZ]; double PHIC2y[NX][NY][NZ];
	double PHIC1z[NX][NY][NZ]; double PHIC2z[NX][NY][NZ];
	//-------------------------------------------------------------------------------------
	


	//-------------------------------------------------------------------------------------
	for ( unsigned int ex=0; ex<2; ex++ ){          // 1
		for ( unsigned int ey=0; ey<2; ey++ ){      // 22
			for ( unsigned int ez=0; ez<2; ez++ ){  // 333
				// For each combination eee, eeo, eoe, eoo, ..., ooo:
				//-------------------------------------------------------------------------
				// x-filtering (pad with zeros where needed)
				v[0] = f[ex][0]; v[1] = f[ex][1]; v[2] = f[ex][2];
				for ( unsigned int m=0; m<NZ; m++ ){  // z-coordinate
					for ( unsigned int l=0; l<NY; l++ ){ // y-coordinate
						PHIC1x[0][l][m] = (v[1])*(PHIx[0][l][m]) + (v[2])*(PHIx[1][l][m]);
						PHIC1y[0][l][m] = (v[1])*(PHIy[0][l][m]) + (v[2])*(PHIy[1][l][m]);
						PHIC1z[0][l][m] = (v[1])*(PHIz[0][l][m]) + (v[2])*(PHIz[1][l][m]);
						for ( unsigned int k=1; k<NX-1; k++ ){ // x-coordinate
							PHIC1x[k][l][m] = (v[0])*(PHIx[k-1][l][m]) + (v[1])*(PHIx[k][l][m]) + (v[2])*(PHIx[k+1][l][m]);
							PHIC1y[k][l][m] = (v[0])*(PHIy[k-1][l][m]) + (v[1])*(PHIy[k][l][m]) + (v[2])*(PHIy[k+1][l][m]);
							PHIC1z[k][l][m] = (v[0])*(PHIz[k-1][l][m]) + (v[1])*(PHIz[k][l][m]) + (v[2])*(PHIz[k+1][l][m]);
						}
						PHIC1x[NX-1][l][m] = (v[0])*(PHIx[NX-2][l][m]) + (v[1])*(PHIx[NX-1][l][m]);
						PHIC1y[NX-1][l][m] = (v[0])*(PHIy[NX-2][l][m]) + (v[1])*(PHIy[NX-1][l][m]);
						PHIC1z[NX-1][l][m] = (v[0])*(PHIz[NX-2][l][m]) + (v[1])*(PHIz[NX-1][l][m]);
					}
				}
				//-------------------------------------------------------------------------
				// y-filtering  (pad with zeros where needed)
				v[0] = f[ey][0]; v[1] = f[ey][1]; v[2] = f[ey][2];
				for ( unsigned int m=0; m<NZ; m++ ){  // z-coordinate
					for ( unsigned int k=0; k<NX; k++ ){ // x-coordinate
						PHIC2x[k][0][m] = (v[1])*(PHIC1x[k][0][m]) + (v[2])*(PHIC1x[k][1][m]);
						PHIC2y[k][0][m] = (v[1])*(PHIC1y[k][0][m]) + (v[2])*(PHIC1y[k][1][m]);
						PHIC2z[k][0][m] = (v[1])*(PHIC1z[k][0][m]) + (v[2])*(PHIC1z[k][1][m]);
						for ( unsigned int l=1; l<NY-1; l++ ){ // y-coordinate
							PHIC2x[k][l][m] = (v[0])*(PHIC1x[k][l-1][m]) + (v[1])*(PHIC1x[k][l][m]) + (v[2])*(PHIC1x[k][l+1][m]);
							PHIC2y[k][l][m] = (v[0])*(PHIC1y[k][l-1][m]) + (v[1])*(PHIC1y[k][l][m]) + (v[2])*(PHIC1y[k][l+1][m]);
							PHIC2z[k][l][m] = (v[0])*(PHIC1z[k][l-1][m]) + (v[1])*(PHIC1z[k][l][m]) + (v[2])*(PHIC1z[k][l+1][m]);
						}
						PHIC2x[k][NY-1][m] = (v[0])*(PHIC1x[k][NY-2][m]) + (v[1])*(PHIC1x[k][NY-1][m]);
						PHIC2y[k][NY-1][m] = (v[0])*(PHIC1y[k][NY-2][m]) + (v[1])*(PHIC1y[k][NY-1][m]);
						PHIC2z[k][NY-1][m] = (v[0])*(PHIC1z[k][NY-2][m]) + (v[1])*(PHIC1z[k][NY-1][m]);
					}
				}
				//-------------------------------------------------------------------------
				// z-filtering  (pad with zeros where needed)
				v[0] = f[ez][0]; v[1] = f[ez][1]; v[2] = f[ez][2];
				for ( unsigned int k=0; k<2*NX-2; k+=2 ){  // x-coordinate
					if( (k==2*NX-4) && (ex==1) ) break;
					for ( unsigned int l=0; l<2*NY-2; l+=2 ){ // y-coordinate
						if( (l==2*NY-4) && (ey==1) ) break;
						for ( unsigned int m=0; m<2*NZ-2; m+=2 ){ // z-coordinate
							if( (m==2*NZ-4) && (ez==1) ) break;
							PHI2x[k+ex][l+ey][m+ez] = (v[0])*(PHIC2x[k/2+ex][l/2+ey][m/2+ez-index[ez]]) +
								(v[1])*(PHIC2x[k/2+ex][l/2+ey][m/2+ez]) + (v[2])*(PHIC2x[k/2+ex][l/2+ey][m/2+ez+1]);
							PHI2y[k+ex][l+ey][m+ez] = (v[0])*(PHIC2y[k/2+ex][l/2+ey][m/2+ez-index[ez]]) +
								(v[1])*(PHIC2y[k/2+ex][l/2+ey][m/2+ez]) + (v[2])*(PHIC2y[k/2+ex][l/2+ey][m/2+ez+1]);
							PHI2z[k+ex][l+ey][m+ez] = (v[0])*(PHIC2z[k/2+ex][l/2+ey][m/2+ez-index[ez]]) +
								(v[1])*(PHIC2z[k/2+ex][l/2+ey][m/2+ez]) + (v[2])*(PHIC2z[k/2+ex][l/2+ey][m/2+ez+1]);
						}
					}
				}
				//-------------------------------------------------------------------------
			} // 333
		} // 22
	} // 1
	//-------------------------------------------------------------------------------------


	
	//-------------------------------------------------------------------------------------
	// Free memory asociated with old buffer
	for(p=0; p<NX; p++){
		for(q=0; q<NY; q++){
			delete PHIx[p][q];
			delete PHIy[p][q];
			delete PHIz[p][q];
		}
		delete PHIx[p];
		delete PHIy[p];
		delete PHIz[p];
	}
	delete PHIx;
	delete PHIy;
	delete PHIz;
	//-------------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------------
	// Update the buffer of spline coefficients:
	PHIx = PHI2x;
	PHIy = PHI2y;
	PHIz = PHI2z;
	//-------------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------------
	// Update the buffer for Jacobian:
	m_Jacobian.set_size(3,3*(2*NX-3)*(2*NY-3)*(2*NZ-3));
	m_Jacobian.Fill(NumericTraits<typename JacobianType::ValueType>::Zero);
	// Update the values of NX and NY to the next level of resolution:
	NX = 2*NX - 3;
	NY = 2*NY - 3;
	NZ = 2*NZ - 3;
	//-------------------------------------------------------------------------------------
	
	return;
}

/** Interpolate given points */ //[20-05-2005-12:40]
template <class TScalarType, unsigned int NDimensions>
void D3BSplineTransform<TScalarType, NDimensions>
::InterpolatePoints(std::vector<InputPointType> points, std::vector<OutputPointType> values)
{
	/* This method assumes you want to interpolate given DISPLACEMENTS 'values' at given
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
	// Counters:
	unsigned int k;
	unsigned int l;
	unsigned int m;
	// Number of points for interpolation:
	unsigned int NP = points.size();
	// The PHI coefficients:
	double PHI[4][4][4];
	// The divider coefficients:
	double OMEGA[NX][NY][NZ];
	// The weighting factors:
	double DELTAx[NX][NY][NZ]; 
	double DELTAy[NX][NY][NZ];
	double DELTAz[NX][NY][NZ];
	// Zero the dividers and weighting factors:
	for(k=0; k<NX; k++){
		for(l=0; l<NY; l++){
			for(m=0; m<NZ; m++){
				OMEGA[k][l][m]  = 0.0;
				DELTAx[k][l][m] = 0.0; 
				DELTAy[k][l][m] = 0.0;
				DELTAz[k][l][m] = 0.0;
			}
		}
	}
	
	// Current Interpolating point and current Interpolating value:
	InputPointType  cpoint;
	OutputPointType cvalue;
	// Grid indices:
	unsigned int i;
	unsigned int j;
	unsigned int h;
	// s and t parameters:
	double x;
	double y;
	double z;
	// wkl's sum:
	double wsq;

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
		z = ((double)(NZ-3))*( ( (double)(cpoint[2]) - m_Origin[2] ) / m_Extent[2] );
		
		// If the point is outside the region the grid is defined on, we ignore that point:
		if( !(x<0 || x>=NX-3 || y<0 || y>=NY-3 || z<0 || z>=NZ-3) )
		{
			//  2.- Get the i and j indices:
			i = (unsigned int)(floor(x));
			j = (unsigned int)(floor(y));
			h = (unsigned int)(floor(z));
			//  3.- Get the s and t parameters:
			x = x - floor(x);
			y = y - floor(y);
			z = z - floor(z);

			// Compute \phi_{kl}: and add it to the corresponding \delta and \omega locations
			for( k=0; k<4; k++ ){
				for( l=0; l<4; l++ ){
					for( m=0; m<4; m++ ){
						PHI[k][l][m] = Bi(x,k)*Bi(y,l)*Bi(z,m);
					}
				}
			}
			// Compute wkl's sum:
			wsq = 0.0;
			for(k=0; k<4; k++){
				for(l=0; l<4; l++){
					for(m=0; m<4; m++)
						wsq += (PHI[k][l][m])*(PHI[k][l][m]);
				}
			}
			// Compute \omega and \delta:
			for(k=0; k<4; k++)
			{
				for(l=0; l<4; l++)
				{
					for(m=0; m<4; m++){
						DELTAx[i+k][j+l][h+m] += (PHI[k][l][m])*(PHI[k][l][m])*(((PHI[k][l][m])*((double)(cvalue[0])))/wsq);
						DELTAy[i+k][j+l][h+m] += (PHI[k][l][m])*(PHI[k][l][m])*(((PHI[k][l][m])*((double)(cvalue[1])))/wsq);
						DELTAz[i+k][j+l][h+m] += (PHI[k][l][m])*(PHI[k][l][m])*(((PHI[k][l][m])*((double)(cvalue[2])))/wsq);
						OMEGA[i+k][j+l][h+m]  += (PHI[k][l][m])*(PHI[k][l][m]);
					}
				}
			}
		}
	}

	// For all i,j,h do:
	for(k=0; k<NX; k++){
		for(l=0; l<NY; l++){
			for(m=0; m<NZ; m++){
				if( OMEGA[k][l][m] > 1e-12 ){
					PHIx[k][l][m] += (DELTAx[k][l][m]) / (OMEGA[k][l][m]);
					PHIy[k][l][m] += (DELTAy[k][l][m]) / (OMEGA[k][l][m]);
					PHIz[k][l][m] += (DELTAz[k][l][m]) / (OMEGA[k][l][m]);
				}
			}
		}
	}

	return;
}









/** Interpolate given displacement image */
template <class TScalarType, unsigned int NDimensions>
void D3BSplineTransform<TScalarType, NDimensions>::InterpolateImage( InterpolatedImageType* image )
{
	// Counters:
	unsigned int k;
	unsigned int l;
	unsigned int m;
	// The PHI coefficients:
	double PHI[4][4][4];
	// The divider coefficients:
	double OMEGA[NX][NY][NZ];
	// The weighting factors:
	double DELTAx[NX][NY][NZ]; 
	double DELTAy[NX][NY][NZ];
	double DELTAz[NX][NY][NZ];
	// Zero the dividers and weighting factors:
	for(k=0; k<NX; k++){
		for(l=0; l<NY; l++){
			for(m=0; m<NZ; m++){
				OMEGA[k][l][m]  = 0.0;
				DELTAx[k][l][m] = 0.0; 
				DELTAy[k][l][m] = 0.0;
				DELTAz[k][l][m] = 0.0;
			}
		}
	}
	
	// Current Interpolating point and current Interpolating value:
	typename InterpolatedImageType::PointType cpoint;
	typename InterpolatedImageType::PixelType cvalue;
	OutputPointType op;
	// Grid indices:
	unsigned int i;
	unsigned int j;
	unsigned int h;
	// s and t parameters:
	double x;
	double y;
	double z;
	// wkl's sum:
	double wsq;

	// Image iterator:
	itk::ImageRegionIteratorWithIndex< InterpolatedImageType > it( image, image->GetRequestedRegion() );
	it.GoToBegin();

	// Covert from index to phisical point:
	typename InterpolatedImageType::PointType origin      = image->GetOrigin();
	SpacingType spacing                                   = image->GetSpacing();

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
		z = ((double)(NZ-3))*( ( (double)(cpoint[2]) - m_Origin[2] ) / m_Extent[2] );
		
		// If the point is outside the region the grid is defined on, we ignore that point:
		if( !(x<0 || x>=NX-3 || y<0 || y>=NY-3 || z<0 || z>=NZ-3) )
		{
			//  2.- Get the i and j indices:
			i = (unsigned int)(floor(x));
			j = (unsigned int)(floor(y));
			h = (unsigned int)(floor(z));
			//  3.- Get the s and t parameters:
			x = x - floor(x);
			y = y - floor(y);
			z = z - floor(z);
			
			// Compute \phi_{klm}: and add it to the corresponding \delta and \omega locations
			for( k=0; k<4; k++ ){
				for( l=0; l<4; l++ ){
					for( m=0; m<4; m++ )
						PHI[k][l][m] = Bi(x,k)*Bi(y,l)*Bi(z,m);
				}
			}

			// Compute wkl's sum:
			wsq = 0.0;
			for(k=0; k<4; k++){
				for(l=0; l<4; l++){
					for(m=0; m<4; m++)
						wsq += (PHI[k][l][m])*(PHI[k][l][m]);
				}
			}
			// Compute \omega and \delta:
			for(k=0; k<4; k++)
			{
				for(l=0; l<4; l++)
				{
					for(m=0; m<4; m++){
						DELTAx[i+k][j+l][h+m] += (PHI[k][l][m])*(PHI[k][l][m])*(((PHI[k][l][m])*((double)(cvalue[0])))/wsq);
						DELTAy[i+k][j+l][h+m] += (PHI[k][l][m])*(PHI[k][l][m])*(((PHI[k][l][m])*((double)(cvalue[1])))/wsq);
						DELTAz[i+k][j+l][h+m] += (PHI[k][l][m])*(PHI[k][l][m])*(((PHI[k][l][m])*((double)(cvalue[2])))/wsq);
						OMEGA[i+k][j+l][h+m]  += (PHI[k][l][m])*(PHI[k][l][m]);
					}
				}
			}
		}
		++it;
	}

	// For all i,j do:
	for(k=0; k<NX; k++){
		for(l=0; l<NY; l++){
			for(m=0; m<NZ; m++){
				if( OMEGA[k][l][m] > 1e-12 ){
					PHIx[k][l][m] -= (DELTAx[k][l][m]) / (OMEGA[k][l][m]);
					PHIy[k][l][m] -= (DELTAy[k][l][m]) / (OMEGA[k][l][m]);
					PHIz[k][l][m] -= (DELTAz[k][l][m]) / (OMEGA[k][l][m]);
				}
			}
		}
	}

	return;
}









































// Get the parameters [20-5-2005-9:45]
template <class TScalarType, unsigned int NDimensions>
const typename D3BSplineTransform<TScalarType, NDimensions>::ParametersType &
D3BSplineTransform<TScalarType, NDimensions>
::GetParameters( void ) const
{
	// Appropriatelly resize parameters buffer:
	(*m_ReturnedParameters).SetSize(3*NX*NY*NZ);
	// x-coordinate spline coefficients:
	for( unsigned int k=0; k<NX; k++ ){
		for( unsigned int l=0; l<NY; l++ ){
			for( unsigned int m=0; m<NZ; m++ )
				(*m_ReturnedParameters)[NX*NY*m + NX*l + k] = PHIx[k][l][m];
		}
	}

	// y-coordinate spline coefficients:
	for( unsigned int k=0; k<NX; k++ ){
		for( unsigned int l=0; l<NY; l++ ){
			for( unsigned int m=0; m<NZ; m++ )
				(*m_ReturnedParameters)[NX*NY*NZ + NX*NY*m + NX*l + k] = PHIy[k][l][m];
		}
	}

	// z-coordinate spline coefficients:
	for( unsigned int k=0; k<NX; k++ ){
		for( unsigned int l=0; l<NY; l++ ){
			for( unsigned int m=0; m<NZ; m++ )
				(*m_ReturnedParameters)[2*NX*NY*NZ + NX*NY*m + NX*l + k] = PHIz[k][l][m];
		}
	}

	return (*m_ReturnedParameters);
}





// Set the parameters [20-5-2005-9:50]
template <class TScalarType, unsigned int NDimensions>
void D3BSplineTransform<TScalarType, NDimensions>
::SetParameters( const ParametersType & parameters )
{
	/* As in D3BSplineDeformableTransform, we keep only a pointer to
	the externally managed parameters, but not the actual parameters, 
	because they are already in the PHI matrices*/
	
	if( parameters.Size() != 3*NX*NY*NZ ){
		itkExceptionMacro(<< "Unexpected number of parameters. Needed: " << 3*NX*NY*NZ << 
		               ", provided: " << parameters.Size() << ". Make sure you have called SetNumberOfPoints().");
	}
	
	m_ReturnedParameters = (itk::Array<double>*)(&parameters);	
	
	for( unsigned int k=0; k<NX; k++ ){
		for( unsigned int l=0; l<NY; l++ ){
			for ( unsigned int m=0; m<NZ; m++ ){
				PHIx[k][l][m] = parameters[NX*NY*m + NX*l + k];              // x-coordinate
				PHIy[k][l][m] = parameters[NX*NY*NZ + NX*NY*m + NX*l + k];   // y-coordinate
				PHIz[k][l][m] = parameters[2*NX*NY*NZ + NX*NY*m + NX*l + k]; // z-coordinate
			}
		}
	}	
	
	return;
}


// Set the parameters so the spline represent an affine transform [23-5-2005-10:45]
template <class TScalarType, unsigned int NDimensions>
void
D3BSplineTransform<TScalarType, NDimensions>
::SetAffineTransform( const ParametersType & parameters )
{
	// In this method, we have to take into account that with 1x1 grid we normalize the input coordinates
	// so they fall into the interval [0,1), what leads to an extra transformation which makes us to modify
	// the base parameters
	
	if( (NX != 4) || (NY != 4) || (NZ !=4) )
		itkExceptionMacro( << "This method is only supported for 1x1x1 grids, and current grid is " << NX-3 << "x" << NY-3 << "x" << NZ-3);
	if(parameters.Size() != 12)
		itkExceptionMacro( << "You must provide 12 parameters for an affine transform");
	
	// Get the affine transform parameters:
	double l11 = parameters[0];
	double l12 = parameters[1];
	double l13 = parameters[2];
	
	double l21 = parameters[3];
	double l22 = parameters[4];
	double l23 = parameters[5];

	double l31 = parameters[6];
	double l32 = parameters[7];
	double l33 = parameters[8];

	// The following combination is necesary because ITK::AffineTransform applies the translation
	// before the rotation and scaling:
	double T1  = parameters[9]  + l11*m_Origin[0] + l12*m_Origin[1] + l13*m_Origin[2] - m_Origin[0];
	double T2  = parameters[10] + l21*m_Origin[0] + l22*m_Origin[1] + l23*m_Origin[2] - m_Origin[1];
	double T3  = parameters[11] + l31*m_Origin[0] + l32*m_Origin[1] + l33*m_Origin[2] - m_Origin[2];
	
	// Get the affine transform parameters:
	l11 = l11*m_Extent[0] - m_Extent[0] + 1.0;
	l12 = l12*m_Extent[1];
	l13 = l13*m_Extent[2];

	l21 = l21*m_Extent[0];
	l22 = l22*m_Extent[1] - m_Extent[1] + 1.0;
	l23 = l23*m_Extent[2];

	l31 = l31*m_Extent[0];
	l32 = l32*m_Extent[1];
	l33 = l33*m_Extent[2] - m_Extent[2] + 1.0;
	
	double PHIxBase = T1 - l11 - l12 - l13 + 1;
	double PHIyBase = T2 - l21 - l22 - l23 + 1;
	double PHIzBase = T3 - l31 - l32 - l33 + 1;
	
	for (unsigned int k=0; k<4; k++){
		for (unsigned int l=0; l<4; l++){
			for (unsigned int m=0; m<4; m++){
				PHIx[k][l][m] = PHIxBase + ((double)k)*(l11 - 1.0) + ((double)l)*l12 + ((double)m)*l13;
				PHIy[k][l][m] = PHIyBase + ((double)l)*(l22 - 1.0) + ((double)k)*l21 + ((double)m)*l23;
				PHIz[k][l][m] = PHIzBase + ((double)m)*(l33 - 1.0) + ((double)k)*l31 + ((double)l)*l32;
			}
		}
	}
	
	return;
}



// Print self [20-5-2005-10:00]
template<class TScalarType, unsigned int NDimensions>
void
D3BSplineTransform<TScalarType, NDimensions>::
PrintSelf(std::ostream &os, Indent indent) const 
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "NX: " << NX << std::endl;
  os << indent << "NY: " << NY << std::endl;
  os << indent << "NZ: " << NZ << std::endl;
  os << indent << "m_Origin = [" << m_Origin[0] << ", " << m_Origin[1] << ", " << m_Origin[2] << "]" << std::endl;
  os << indent << "m_Extent = [" << m_Extent[0] << ", " << m_Extent[1] << ", " << m_Extent[2] << "]" << std::endl;
  
  os << indent << "PHIx = " << std::endl;
  for(unsigned int m=0; m<NZ; m++){
	  os << indent << "(Dimension z=" << m << ")" << std::endl;
	  os << indent << "[" << std::endl;
	  for(unsigned int k=0; k<NX; k++){
		  os << indent;
		  for(unsigned int l=0; l<NY; l++)
			  os << PHIx[k][l][m] << ", ";
		  os << std::endl;
	  }
	  os << indent << "]" << std::endl;
  }

  os << std::endl;

  os << indent << "PHIy = " << std::endl;
  for(unsigned int m=0; m<NZ; m++){
	  os << indent << "(Dimension z=" << m << ")" << std::endl;
	  os << indent << "[" << std::endl;
	  for(unsigned int k=0; k<NX; k++){
		  os << indent;
		  for(unsigned int l=0; l<NY; l++)
			  os << PHIy[k][l][m] << ", ";
		  os << std::endl;
	  }
	  os << indent << "]" << std::endl;
  }

  os << std::endl;

  os << indent << "PHIz = " << std::endl;
  for(unsigned int m=0; m<NZ; m++){
	  os << indent << "(Dimension z=" << m << ")" << std::endl;
	  os << indent << "[" << std::endl;
	  for(unsigned int k=0; k<NX; k++){
		  os << indent;
		  for(unsigned int l=0; l<NY; l++)
			  os << PHIz[k][l][m] << ", ";
		  os << std::endl;
	  }
	  os << indent << "]" << std::endl;
  }
}


// Transform a point [20-5-2005-12:50]
template<class TScalarType, unsigned int NDimensions>
typename D3BSplineTransform<TScalarType, NDimensions>::OutputPointType
D3BSplineTransform<TScalarType, NDimensions>::
TransformPoint(const InputPointType &point) const 
{
	// Offset in both directions:
	InputVectorType offset;
	// Displacement along each coordinate:
	double xoffset = 0.0;
	double yoffset = 0.0;
	double zoffset = 0.0;
	// Grid indices:
	unsigned int i;
	unsigned int j;
	unsigned int h;
	// Conversion to double type:
	double x = (double)(point[0]);
	double y = (double)(point[1]);	
	double z = (double)(point[2]);

	// If the point is outside the region the grid is defined on, we return the same point:
	if( x<m_Origin[0] || x>m_Origin[0]+m_Extent[0] || y<m_Origin[1] || y>m_Origin[1]+m_Extent[1] || z<m_Origin[2] || z>m_Origin[2]+m_Extent[2])
		return point;

	// Normalization of coordinates to grid coordinates:
	x = ((double)(NX-3))*( ( x - m_Origin[0] ) / m_Extent[0] );
	y = ((double)(NY-3))*( ( y - m_Origin[1] ) / m_Extent[1] );
	z = ((double)(NZ-3))*( ( z - m_Origin[2] ) / m_Extent[2] );
	
	// Get the i and j indices:
	i = (unsigned int)(floor(x));
	j = (unsigned int)(floor(y));
	h = (unsigned int)(floor(z));
	
	// Get the s and t parameters:
	x = x - floor(x);
	y = y - floor(y);
	z = z - floor(z);

	// Get the displacements by means of the evaluation of splines:
	for( unsigned int k=0; k<4; k++ ){
		for( unsigned int l=0; l<4; l++ ){
			for( unsigned int m=0; m<4; m++ ){
				xoffset += (PHIx[i+k][j+l][h+m])*Bi(x,k)*Bi(y,l)*Bi(z,m);
				yoffset += (PHIy[i+k][j+l][h+m])*Bi(x,k)*Bi(y,l)*Bi(z,m);
				zoffset += (PHIz[i+k][j+l][h+m])*Bi(x,k)*Bi(y,l)*Bi(z,m);
			}
		}
	}

	// Update the value of the point:
	offset[0] = xoffset;
	offset[1] = yoffset;
	offset[2] = zoffset;
	
	// return:
	return (point + offset);
}
  




// Compute the Jacobian in one position [20-5-2005-13:05]
template<class TScalarType, unsigned int NDimensions>
const typename D3BSplineTransform<TScalarType, NDimensions>::JacobianType & 
D3BSplineTransform< TScalarType, NDimensions >::
GetJacobian( const InputPointType &point) const
{
	unsigned int k;
	unsigned int l;
	unsigned int m;

	unsigned int i;
	unsigned int j;
	unsigned int h;

	// Conversion to double type:
	double x = (double)(point[0]);
	double y = (double)(point[1]);
	double z = (double)(point[2]);
	
	// First, we set to 0 the values of Jacobian, when needed:
	for(k=m_LastVisited[0]; k<m_LastVisited[0]+4; k++){
		for(l=m_LastVisited[1]; l<m_LastVisited[1]+4; l++){
			for(m=m_LastVisited[2]; m<m_LastVisited[2]+4; m++){
				m_Jacobian(0,k + l*NX + m*NX*NY)              = 0.0;
				m_Jacobian(1,k + l*NX + m*NX*NY + NX*NY*NZ)   = 0.0;
				m_Jacobian(2,k + l*NX + m*NX*NY + 2*NX*NY*NZ) = 0.0;
			}
		}
	}

	// If the point is outside the region the grid is defined on, we return 0:
	if( x<m_Origin[0] || x>m_Origin[0]+m_Extent[0] || y<m_Origin[1] || y>m_Origin[1]+m_Extent[1] || z<m_Origin[2] || z>m_Origin[2]+m_Extent[2])
		return m_Jacobian;

	// Normalization of coordinates to grid coordinates:
	x = ((double)(NX-3))*( ( x - m_Origin[0] ) / m_Extent[0] );
	y = ((double)(NY-3))*( ( y - m_Origin[1] ) / m_Extent[1] );
	z = ((double)(NZ-3))*( ( z - m_Origin[2] ) / m_Extent[2] );
	
	// Get the i and j indices:
	i = (unsigned int)(floor(x));
	j = (unsigned int)(floor(y));
	h = (unsigned int)(floor(z));
	
	// Get the s and t parameters:
	x = x - floor(x);
	y = y - floor(y);
	z = z - floor(z);

	// Update the LastVisited:
	m_LastVisited[0] = i;
	m_LastVisited[1] = j;
	m_LastVisited[2] = h;
	
	// Fill the rest of Jacobian values:
	for( k=i; k<i+4; k++ ){
		for( l=j; l<j+4; l++ ){
			for( m=h; m<h+4; m++ ){
				m_Jacobian(0,k + l*NX + m*NX*NY)              = Bi(x,k-i)*Bi(y,l-j)*Bi(z,m-h);
				m_Jacobian(1,k + l*NX + m*NX*NY + NX*NY*NZ)   = Bi(x,k-i)*Bi(y,l-j)*Bi(z,m-h);
				m_Jacobian(2,k + l*NX + m*NX*NY + 2*NX*NY*NZ) = Bi(x,k-i)*Bi(y,l-j)*Bi(z,m-h);
			}
		}
	}

	// return:
	return m_Jacobian;
}


// Compute the penalty term
template<class TScalarType, unsigned int NDimensions>
double D3BSplineTransform< TScalarType, NDimensions >::GetCorrection()
{
	// Initialization ckecking:
	if( !(NX*NY*NZ) )
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
	double nz  = (double)NZ-3.0;

	double X   = (double)(m_Extent[0]);
	double Y   = (double)(m_Extent[1]);
	double Z   = (double)(m_Extent[2]);
	
	double la1 = (nx*nx*nx)/(ny*nz*X*X*X*X);
	double la2 = (ny*ny*ny)/(nx*nz*Y*Y*Y*Y);
	double la3 = (nz*nz*nz)/(nx*ny*Z*Z*Z*Z);
	
	double la4 = 2.0*(nx*ny)/(nz*X*X*Y*Y);
	double la5 = 2.0*(ny*nz)/(nx*Y*Y*Z*Z);
	double la6 = 2.0*(nx*nz)/(ny*X*X*Z*Z);
	
	// Reckoning:
	double val = 0.0;
	for( unsigned int i=0; i<NX-3; i++ ){ // 1
		for ( unsigned int j=0; j<NY-3; j++ ){ // 22
			for ( unsigned int h=0; h<NZ-3; h++ ){ // 333
				// For each parameter (i, j), compute the integral in the cell:
				for( unsigned int k=0; k<4; k++ ){ // 4444
					for( unsigned int l=0; l<4; l++ ){ // 55555
						for( unsigned int m=0; m<4; m++ ){ // 666666
						// First sums
							for( unsigned int p=0; p<4; p++ ){ // 7777777
								for( unsigned int q=0; q<4; q++ ){ // 88888888
									for( unsigned int r=0; r<4; r++ ){ // 999999999
										// Second sums
										val += ( (PHIx[i+k][j+l][h+m])*(PHIx[i+p][j+q][h+r]) + (PHIy[i+k][j+l][h+m])*(PHIy[i+p][j+q][h+r])
											+ (PHIy[i+k][j+l][h+m])*(PHIy[i+p][j+q][h+r]) )*
											( la1*(D2D2[k][p])*(D0D0[l][q])*(D0D0[m][r]) + la2*(D0D0[k][p])*(D2D2[l][q])*(D0D0[m][r]) +
											  la3*(D0D0[k][p])*(D0D0[l][q])*(D2D2[m][r]) + la4*(D1D1[k][p])*(D1D1[l][q])*(D0D0[m][r]) +
											  la5*(D0D0[k][p])*(D1D1[l][q])*(D1D1[m][r]) + la6*(D1D1[k][p])*(D0D0[l][q])*(D1D1[m][r]) );
									} // 999999999
								} // 88888888
							} // 7777777
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
typename D3BSplineTransform<TScalarType, NDimensions>::CorrectionJacobianType 
D3BSplineTransform< TScalarType, NDimensions >::GetCorrectionJacobian()
{
	// Initialization ckecking:
	if( !(NX*NY*NZ) )
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
	double nz  = (double)NZ-3.0;

	double X   = (double)(m_Extent[0]);
	double Y   = (double)(m_Extent[1]);
	double Z   = (double)(m_Extent[2]);
	
	double la1 = (nx*nx*nx)/(ny*nz*X*X*X*X);
	double la2 = (ny*ny*ny)/(nx*nz*Y*Y*Y*Y);
	double la3 = (nz*nz*nz)/(nx*ny*Z*Z*Z*Z);
	
	double la4 = 2.0*(nx*ny)/(nz*X*X*Y*Y);
	double la5 = 2.0*(ny*nz)/(nx*Y*Y*Z*Z);
	double la6 = 2.0*(nx*nz)/(ny*X*X*Z*Z);

	// Declare the Jacobian and set the appropriate size:
	CorrectionJacobianType cj; 
	cj.set_size(3*NX*NY*NZ);

	// Cumulative values of each component of Jacobian
	double valx;
	double valy;
	double valz;
	// Define the lower limit of the first sums:
	unsigned int k0;
	unsigned int l0;
	unsigned int m0;
	// Define the upper limit of the first sums:
	unsigned int k3;
	unsigned int l3;
	unsigned int m3;

	// Reckoning of each Jacobian value:
	k0 = 0; k3 = 0;                  // Initiallize k0 and k3
	for( unsigned int r=0; r<NX; r++ ){ // 1
		if( r>=NX-3 ) k0++;          // Update the lower limit of first sumations, if needed
		l0 = 0; l3 = 0;              // Initiallize l0 and l3
		for ( unsigned int s=0; s<NY; s++ ){ // 22
			if( s>=NY-3 ) l0++;      // Update the lower limit of first sumations, if needed
			m0 = 0; m3 = 0;          // Initiallize m0 and m3
			for ( unsigned int t=0; t<NZ; t++){ // 333
				if( t>=NZ-3 ) m0++;  // Update the lower limit of first sumations, if needed
				// For each parameter (r, s, t), compute the Jacobian component by means of a double sumation:
				valx = 0.0;
				valy = 0.0;
				valz = 0.0;
				for( unsigned int k=k0; k<=k3; k++ ){ // 4444
					for( unsigned int l=l0; l<=l3; l++ ){ // 55555
						for( unsigned int m=m0; m<=m3; m++ ){ // 666666
							// sumation over each influence cell
							for( unsigned int p=0; p<=3; p++ ){ // 7777777
								for( unsigned int q=0; q<=3; q++ ){ // 88888888
									for( unsigned int o=0; o<=3; o++ ){ // 999999999
									// Pre-computed integral in each influence cell
									valx += 2.0*(PHIx[p+r-k][q+s-l][o+t-m])*
										( la1*(D2D2[k][p])*(D0D0[l][q])*(D0D0[m][o]) + la2*(D0D0[k][p])*(D2D2[l][q])*(D0D0[m][o]) 
										+ la3*(D0D0[k][p])*(D0D0[l][q])*(D2D2[m][o]) + la4*(D1D1[k][p])*(D1D1[l][q])*(D0D0[m][o])
										+ la5*(D0D0[k][p])*(D1D1[l][q])*(D1D1[m][o]) + la6*(D1D1[k][p])*(D0D0[l][q])*(D1D1[m][o]) );
									valy += 2.0*(PHIy[p+r-k][q+s-l][o+t-m])*
										( la1*(D2D2[k][p])*(D0D0[l][q])*(D0D0[m][o]) + la2*(D0D0[k][p])*(D2D2[l][q])*(D0D0[m][o]) 
										+ la3*(D0D0[k][p])*(D0D0[l][q])*(D2D2[m][o]) + la4*(D1D1[k][p])*(D1D1[l][q])*(D0D0[m][o])
										+ la5*(D0D0[k][p])*(D1D1[l][q])*(D1D1[m][o]) + la6*(D1D1[k][p])*(D0D0[l][q])*(D1D1[m][o]) );
									valz += 2.0*(PHIz[p+r-k][q+s-l][o+t-m])*
										( la1*(D2D2[k][p])*(D0D0[l][q])*(D0D0[m][o]) + la2*(D0D0[k][p])*(D2D2[l][q])*(D0D0[m][o]) 
										+ la3*(D0D0[k][p])*(D0D0[l][q])*(D2D2[m][o]) + la4*(D1D1[k][p])*(D1D1[l][q])*(D0D0[m][o])
										+ la5*(D0D0[k][p])*(D1D1[l][q])*(D1D1[m][o]) + la6*(D1D1[k][p])*(D0D0[l][q])*(D1D1[m][o]) );				
									} // 999999999
								} // 88888888
							} // 7777777
						} // 666666
					} // 55555
				} // 4444
				//-----------------------------------------------------------------------------------------------
				// Once here, we've got already the value of all Jacobian values:
				cj[NX*NY*t+NX*s+r]            = valx;
				cj[NX*NY*NZ+NX*NY*t+NX*s+r]   = valy;
				cj[2*NX*NY*NZ+NX*NY*t+NX*s+r] = valy;
				//-----------------------------------------------------------------------------------------------
				if( t<3 ) m3++; // Update the upper limit of first sumations, if needed
			} // 333		
			if( s<3 ) l3++;     // Update the upper limit of first sumations, if needed
		} // 22
		if( r<3 ) k3++;         // Update the upper limit of first sumations, if needed
	} // 1

	
	// return:
	return cj;
}



// Get the total number of parameters needed [20-5-2005-10:00]
template<class TScalarType, unsigned int NDimensions>
unsigned int D3BSplineTransform<TScalarType, NDimensions>::GetNumberOfParameters(void) const
{
	return 3*NX*NY*NY;
}




// Set the parameters for an Identity transform of this class [20-5-2005-10:00]
template<class TScalarType, unsigned int NDimensions>
void D3BSplineTransform<TScalarType, NDimensions>::
SetIdentity()
{
	if(NX*NY*NZ){
		for(unsigned int k=0; k<NX; k++){
			for(unsigned int l=0; l<NY; l++){
				for(unsigned int m=0; m<NZ; m++){
					PHIx[k][l][m] = 0.0;
					PHIy[k][l][m] = 0.0;
					PHIz[k][l][m] = 0.0;
				}
			}
		}
	}
	return;
}
 
  
} // namespace

#endif
