/*=========================================================================

	Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSaturnTensorGlyph.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSaturnTensorGlyph.h"

#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkSuperquadricSource.h"

// Construct object with scaling on and scale factor 1.0. Eigenvalues are 
// extracted, glyphs are colored with input scalar data, and logarithmic
// scaling is turned off.
vtkSaturnTensorGlyph::vtkSaturnTensorGlyph()
{

	this->input = NULL;
	this->inputPoints = NULL;

	this->Scaling = 1;
	this->ScaleFactor = 1.0;

	this->GlyphType = ELLIPSOID;
	this->ColorMode = COLOR_BY_FA;
	this->PhiResolution = 8;
	this->ThetaResolution = 8;

	this->Bounds[0] = 0;
	this->Bounds[1] = 0;
	this->Bounds[2] = 0;
	this->Bounds[3] = 0;
	this->Bounds[4] = 0;
	this->Bounds[5] = 0;

	this->Gamma = 3.0;

	this->FilterMode = FILTER_BY_FA;
	this->FilterThreshold = 0.0;
}

vtkSaturnTensorGlyph::~vtkSaturnTensorGlyph()
{
}

vtkPolyData *vtkSaturnTensorGlyph::GetOutput()
{

	vtkPolyData *output = vtkPolyData::New();

	vtkDataArray *inTensors;
	double tensor[9];
	vtkDataArray *inScalars;
	vtkIdType numPts, numSourcePts, numSourceCells, inPtId;
	vtkPoints *sourcePts;
	vtkDataArray *sourceNormals;
	vtkCellArray *sourceCells, *cells;
	vtkPoints *newPts;
	vtkFloatArray *newScalars=NULL;
	vtkFloatArray *newNormals=NULL;
	double x[3], s;
	vtkTransform *trans;
	vtkCell *cell;
	vtkIdList *cellPts;
	int npts;
	vtkIdType *pts;
	vtkIdType ptIncr, cellId;
	vtkIdType subIncr;
	int numDirs, dir, eigen_dir, symmetric_dir;
	vtkMatrix4x4 *matrix;
	float *m[3], w[3], *v[3];
	float m0[3], m1[3], m2[3];
	float v0[3], v1[3], v2[3];
	float xv[3], yv[3], zv[3];
	float maxScale;
	vtkPointData *pd, *outPD;

	TensorImageType::IndexType pixelIndex;
	TensorPixelType pixel, pixel1, pixel2, pixel3, pixel4, pixel5, pixel6, pixel7, pixel8;
	const TensorPixelType pixel_prueba;
	EigenValuesArrayType eigval;
	EigenVectorsMatrixType eigvec;
	double factor, norm_factor;
	float scalarValue;
	int i,j,k,p;
	int xSize, ySize, zSize;
	RealType cl, cp, cs, t1, t2, t3, num;
	int cont=0;
	double alpha, beta;

	vtkSuperquadricSource *superquadric;

	TensorImageType::PointType origin = this->input->GetOrigin();
	TensorImageType::SpacingType spacing = this->input->GetSpacing();

	switch (this->GlyphType)
	{
		case ELLIPSOID: { 
			vtkSphereSource *sphere = vtkSphereSource::New();
			sphere->SetPhiResolution(this->PhiResolution);
			sphere->SetThetaResolution(this->ThetaResolution);
			this->source = sphere->GetOutput();
			this->source->Update();
//			sphere->Delete();
			break;
		}

		case CUBOID: {
			vtkCubeSource *cube = vtkCubeSource::New();
			source = cube->GetOutput();
			this->source->Update();
			break;
		}

		case SUPERQUADRIC: {
			superquadric = vtkSuperquadricSource::New();
			superquadric->SetPhiResolution(this->PhiResolution);
			superquadric->SetThetaResolution(this->ThetaResolution);
			this->source = superquadric->GetOutput();
			this->source->Update();
			break;
		}
	}

	pts = new vtkIdType[source->GetMaxCellSize()];
	trans = vtkTransform::New();
	matrix = vtkMatrix4x4::New();

	// set up working matrices
	m[0] = m0; m[1] = m1; m[2] = m2; 
	v[0] = v0; v[1] = v1; v[2] = v2; 

	outPD = output->GetPointData();

	xSize = Bounds[1]-Bounds[0] + 1;
	ySize = Bounds[3]-Bounds[2] + 1;
	zSize = Bounds[5]-Bounds[4] + 1;
	numPts = xSize * ySize * zSize;

	//
	// Allocate storage for output PolyData
	//
	sourcePts = source->GetPoints();
	numSourcePts = sourcePts->GetNumberOfPoints();
	numSourceCells = source->GetNumberOfCells();

	newPts = vtkPoints::New();

	// only copy scalar data through
	pd = source->GetPointData();

	newScalars = vtkFloatArray::New();

	if ( (sourceNormals = pd->GetNormals()) )
	{
		newNormals = vtkFloatArray::New();
		newNormals->SetNumberOfComponents(3);
	}


	//
	// Traverse all Input points, transforming glyph at Source points
	//
	trans->PreMultiply();

	for (i=this->Bounds[0]; i<=this->Bounds[1]; i++)
	for (j=this->Bounds[2]; j<=this->Bounds[3]; j++)
	for (k=this->Bounds[4]; k<=this->Bounds[5]; k++)
	{
		ptIncr = inPtId * numSourcePts;

		pixelIndex[0] = i;
		pixelIndex[1] = j;
		pixelIndex[2] = k;
		pixel = input->GetPixel(pixelIndex);

		pixel.ComputeEigenSystem(eigval,eigvec);

		if (eigval[0]==0) continue;

//		if ( (this->GlyphType==SUPERQUADRIC) || (this->ColorMode==COLOR_BY_CL) )
		pixel.ComputeShapeCoefficients(cl,cp,cs);

		switch (this->FilterMode)
		{
			case FILTER_BY_FA: 
				if (pixel.GetFractionalAnisotropy() < FilterThreshold) continue;
				break;

			case FILTER_BY_CL: 
				if (cl < FilterThreshold) continue;
				break;

			case FILTER_BY_CS: 
				if (cs > FilterThreshold) continue;
				break;

		}

		// Now do the real work for each "direction"

		// Remove previous scales ...
		trans->Identity();

		// translate Source to Input point
		x[0] = origin[0] + spacing[0] * pixelIndex[0];
		x[1] = origin[1] + spacing[1] * pixelIndex[1];
		x[2] = origin[2] + spacing[2] * pixelIndex[2];
		trans->Translate(x[0], x[1], x[2]);


		// normalized eigenvectors rotate object for eigen direction 0
		float signo = 1.0;
		if ( (eigvec[0][1]*eigvec[1][2]-eigvec[0][2]*eigvec[1][1]) * eigvec[2][0] < 0) signo = -1.0; 

		matrix->Element[0][0] = eigvec[0][0];
		matrix->Element[0][1] = eigvec[1][0];
		matrix->Element[0][2] = eigvec[2][0];
		matrix->Element[1][0] = eigvec[0][1];
		matrix->Element[1][1] = eigvec[1][1];
		matrix->Element[1][2] = eigvec[2][1];
		matrix->Element[2][0] = signo*eigvec[0][2];
		matrix->Element[2][1] = signo*eigvec[1][2];
		matrix->Element[2][2] = signo*eigvec[2][2];

		trans->Concatenate(matrix);


		if (fabs(eigval[0])>fabs(eigval[2]))
			norm_factor = 1 / eigval[0];

		else norm_factor = 1 / eigval[2];

		factor = norm_factor;

		if (this->Scaling)
			factor *= this->ScaleFactor;

		trans->Scale(factor*eigval[0],factor*eigval[1],factor*eigval[2]);

 
		if (this->GlyphType==SUPERQUADRIC) 
		{
			if (cl<=cp) 
			{
				alpha = pow(1-cp,Gamma);
				beta = pow(1-cl,Gamma);
			}
			else 
			{
				alpha = pow(1-cl,Gamma);
				beta = pow(1-cp,Gamma);
			}

			superquadric->SetThetaRoundness(alpha);
			superquadric->SetPhiRoundness(beta);

			this->source = superquadric->GetOutput();
			this->source->Update();

			sourcePts = source->GetPoints();
			sourceNormals = source->GetPointData()->GetNormals();
	 	}

		// multiply points (and normals if available) by resulting
		// matrix
		trans->TransformPoints(sourcePts,newPts);

		// Apply the transformation to a series of points, 
		// and append the results to outPts.
		if ( newNormals )
		{
			trans->TransformNormals(sourceNormals,newNormals);
		}

		switch (this->ColorMode)
		{
			case COLOR_BY_FA: 
				scalarValue = pixel.GetFractionalAnisotropy();
				break;

			case COLOR_BY_RA: 
				scalarValue = pixel.GetRelativeAnisotropy();
				break;

			case COLOR_BY_CL: 
				scalarValue = cl;
				break;

		}

		// Copy point data from source
		for (p=0; p < numSourcePts; p++) 
		{
			newScalars->InsertNextValue(scalarValue);
		}

		cont++;
	}

	if (this->inputPoints) numPts = this->inputPoints->GetNumberOfPoints();

	if (this->inputPoints)

	for (inPtId=0; inPtId < numPts; inPtId++)
	{
		ptIncr = inPtId * numSourcePts;

		this->inputPoints->GetPoint(inPtId,x);
		interpolacionLogEuclidea(x,&pixel);

		pixel.ComputeEigenSystem(eigval,eigvec);

		if (eigval[0]==0) continue;

		pixel.ComputeShapeCoefficients(cl,cp,cs);

		// Now do the real work for each "direction"

		// Remove previous scales ...
		trans->Identity();

		// translate Source to Input point
		trans->Translate(x[0], x[1], x[2]);


		float signo = 1.0;
		if ( (eigvec[0][1]*eigvec[1][2]-eigvec[0][2]*eigvec[1][1]) * eigvec[2][0] < 0) signo = -1.0; 

		// normalized eigenvectors rotate object for eigen direction 0
		matrix->Element[0][0] = eigvec[0][0];
		matrix->Element[0][1] = eigvec[1][0];
		matrix->Element[0][2] = eigvec[2][0];
		matrix->Element[1][0] = eigvec[0][1];
		matrix->Element[1][1] = eigvec[1][1];
		matrix->Element[1][2] = eigvec[2][1];
		matrix->Element[2][0] = eigvec[0][2];
		matrix->Element[2][1] = eigvec[1][2];
		matrix->Element[2][2] = eigvec[2][2];

		trans->Concatenate(matrix);


		if (fabs(eigval[0])>fabs(eigval[2]))
			norm_factor = 1 / eigval[0];

		else norm_factor = 1 / eigval[2];

		factor = norm_factor;

		if (this->Scaling)
			factor *= this->ScaleFactor;

		trans->Scale(factor*eigval[0],factor*eigval[1],factor*eigval[2]);


		if (this->GlyphType==SUPERQUADRIC) 
		{
			if (cl<=cp) 
			{
				alpha = pow(1-cp,Gamma);
				beta = pow(1-cl,Gamma);
			}
			else 
			{
				alpha = pow(1-cl,Gamma);
				beta = pow(1-cp,Gamma);
			}

			superquadric->SetThetaRoundness(alpha);
			superquadric->SetPhiRoundness(beta);

			this->source = superquadric->GetOutput();
			this->source->Update();

			sourcePts = source->GetPoints();
			sourceNormals = source->GetPointData()->GetNormals();
	 	}

		// multiply points (and normals if available) by resulting
		// matrix
		trans->TransformPoints(sourcePts,newPts);
		// Apply the transformation to a series of points, 
		// and append the results to outPts.
		if ( newNormals )
		{
			trans->TransformNormals(sourceNormals,newNormals);
		}

		switch (this->ColorMode)
		{

			case COLOR_BY_FA: 
				scalarValue = pixel.GetFractionalAnisotropy();
				break;

			case COLOR_BY_RA: 
				scalarValue = pixel.GetRelativeAnisotropy();
				break;

			case COLOR_BY_CL: 
				scalarValue = cl;
				break;

		}

			// Copy point data from source
		for (i=0; i < numSourcePts; i++) 
		{
			newScalars->InsertNextValue(scalarValue);
		}

		cont++;

	}

	// Setting up for calls to PolyData::InsertNextCell()
	if ( (sourceCells=source->GetVerts())->GetNumberOfCells() > 0 )
	{
		cells = vtkCellArray::New();
		cells->Allocate(cont*sourceCells->GetSize());
		output->SetVerts(cells);
		cells->Delete();
	}
	if ( (sourceCells=source->GetLines())->GetNumberOfCells() > 0 )
	{
		cells = vtkCellArray::New();
		cells->Allocate(cont*sourceCells->GetSize());
		output->SetLines(cells);
		cells->Delete();
	}
	if ( (sourceCells=source->GetPolys())->GetNumberOfCells() > 0 )
	{
		cells = vtkCellArray::New();
		cells->Allocate(cont*sourceCells->GetSize());
		output->SetPolys(cells);
		cells->Delete();
	}
	if ( (sourceCells=source->GetStrips())->GetNumberOfCells() > 0 )
	{
		cells = vtkCellArray::New();
		cells->Allocate(cont*sourceCells->GetSize());
		output->SetStrips(cells);
		cells->Delete();
	}

	//
	// First copy all topology (transformation independent)
	//
	for (inPtId=0; inPtId < cont; inPtId++)
	{
		ptIncr = inPtId * numSourcePts;
		for (cellId=0; cellId < numSourceCells; cellId++)
		{
			cell = source->GetCell(cellId);
			cellPts = cell->GetPointIds();
			npts = cellPts->GetNumberOfIds();
			// This variable may be removed, but that 
			// will not improve readability
			subIncr = ptIncr + numSourcePts;
			for (i=0; i < npts; i++)
			{
				pts[i] = cellPts->GetId(i) + subIncr;
			}
			output->InsertNextCell(cell->GetCellType(),npts,pts);
		}
	}

	//
	// Update output and release memory
	//
	delete [] pts;

	output->SetPoints(newPts);
	newPts->Delete();

	if ( newScalars )
	{
		output->GetPointData()->SetScalars(newScalars);
		newScalars->Delete();
	}

	if ( newNormals )
	{
		outPD->SetNormals(newNormals);
		newNormals->Delete();
	}


	trans->Delete();
	matrix->Delete();
//	source->Delete();
	return output;
}

void vtkSaturnTensorGlyph::interpolacionLineal (double x[3],TensorPixelType *pixel) {

		TensorPixelType *puntero;
		TensorPixelType pixel1, pixel2, pixel3, pixel4, pixel5, pixel6, pixel7, pixel8;
		TensorImageType::IndexType pixelIndex;
		TensorImageType::PointType origin = this->input->GetOrigin();
		TensorImageType::SpacingType spacing = this->input->GetSpacing();
		RealType t1,t2,t3;
		
		pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
		pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);
		pixelIndex[2] = (int) ((x[2]-origin[2]) / spacing[2]);

		t1 = ( (x[0]-origin[0]) / spacing[0] - pixelIndex[0] ) / spacing[0];
		t2 = ( (x[1]-origin[1]) / spacing[1] - pixelIndex[1] ) / spacing[1];
		t3 = ( (x[2]-origin[2]) / spacing[2] - pixelIndex[2] ) / spacing[2];

		pixel1 = (input->GetPixel(pixelIndex));

		pixelIndex[0] = pixelIndex[0] + 1;
		pixel2 = (input->GetPixel(pixelIndex));

		pixelIndex[0] = pixelIndex[0] - 1;
		pixelIndex[1] = pixelIndex[1] + 1;
		pixel3 = (input->GetPixel(pixelIndex));

		pixelIndex[0] = pixelIndex[0] + 1;
		pixel4 = (input->GetPixel(pixelIndex));

		pixelIndex[0] = pixelIndex[0] - 1;
		pixelIndex[1] = pixelIndex[1] - 1;
		pixelIndex[2] = pixelIndex[2] + 1;
		pixel5 = (input->GetPixel(pixelIndex));

		pixelIndex[0] = pixelIndex[0] + 1;
		pixel6 = (input->GetPixel(pixelIndex));

		pixelIndex[0] = pixelIndex[0] - 1;
		pixelIndex[1] = pixelIndex[1] + 1;
		pixel7 = (input->GetPixel(pixelIndex));

		pixelIndex[0] = pixelIndex[0] + 1;
		pixel8 = (input->GetPixel(pixelIndex));

		*pixel = ((pixel1*(1-t1) + pixel2*t1) * (1-t2) + (pixel3*(1-t1) + pixel4*t1) * t2) * (1-t3) + ((pixel5*(1-t1) + pixel6*t1) * (1-t2) + (pixel7*(1-t1) + pixel8*t1) * t2) * t3;
	 
}

void vtkSaturnTensorGlyph::interpolacionLogEuclidea (double x[3],TensorPixelType *pixel) {

		TensorPixelType puntero;
		TensorPixelType pixel1, pixel2, pixel3, pixel4, pixel5, pixel6, pixel7, pixel8;
		TensorImageType::IndexType pixelIndex;
		TensorImageType::PointType origin = this->input->GetOrigin();
		TensorImageType::SpacingType spacing = this->input->GetSpacing();
		RealType t1,t2,t3;

		pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
		pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);
		pixelIndex[2] = (int) ((x[2]-origin[2]) / spacing[2]);

		t1 = ( (x[0]-origin[0]) / spacing[0] - pixelIndex[0] ) / spacing[0];
		t2 = ( (x[1]-origin[1]) / spacing[1] - pixelIndex[1] ) / spacing[1];
		t3 = ( (x[2]-origin[2]) / spacing[2] - pixelIndex[2] ) / spacing[2];

		pixel1 = (input->GetPixel(pixelIndex));
		pixel1 = pixel1.ComputeTensorFromLogVector();

		pixelIndex[0] = pixelIndex[0] + 1;
		pixel2 = (input->GetPixel(pixelIndex));
		pixel2 = pixel2.ComputeTensorFromLogVector();

		pixelIndex[0] = pixelIndex[0] - 1;
		pixelIndex[1] = pixelIndex[1] + 1;
		pixel3 = (input->GetPixel(pixelIndex));
		pixel3 = pixel3.ComputeTensorFromLogVector();

		pixelIndex[0] = pixelIndex[0] + 1;
		pixel4 = (input->GetPixel(pixelIndex));
		pixel4 = pixel4.ComputeTensorFromLogVector();

		pixelIndex[0] = pixelIndex[0] - 1;
		pixelIndex[1] = pixelIndex[1] - 1;
		pixelIndex[2] = pixelIndex[2] + 1;
		pixel5 = (input->GetPixel(pixelIndex));
		pixel5 = pixel5.ComputeTensorFromLogVector();

		pixelIndex[0] = pixelIndex[0] + 1;
		pixel6 = (input->GetPixel(pixelIndex));
		pixel6 = pixel6.ComputeTensorFromLogVector();

		pixelIndex[0] = pixelIndex[0] - 1;
		pixelIndex[1] = pixelIndex[1] + 1;
		pixel7 = (input->GetPixel(pixelIndex));
		pixel7 = pixel7.ComputeTensorFromLogVector();

		pixelIndex[0] = pixelIndex[0] + 1;
		pixel8 = (input->GetPixel(pixelIndex));
		pixel8 = pixel8.ComputeTensorFromLogVector();

		*pixel = ((pixel1*(1-t1) + pixel2*t1) * (1-t2) + (pixel3*(1-t1) + pixel4*t1) * t2) * (1-t3) + ((pixel5*(1-t1) + pixel6*t1) * (1-t2) + (pixel7*(1-t1) + pixel8*t1) * t2) * t3;

		*pixel = pixel->ComputeLogVector();

}










