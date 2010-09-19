/*=========================================================================

	Program:   Visualization Toolkit
  Module:    $RCSfile: vtkTensorGlyphDTI.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTensorGlyphDTI.h"

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

// Por defecto, los glifos son elipsoidales con resolución 8, y el color se obtiene de la FA
vtkTensorGlyphDTI::vtkTensorGlyphDTI()
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

vtkTensorGlyphDTI::~vtkTensorGlyphDTI()
{
}

void vtkTensorGlyphDTI::GetOutput(vtkPolyData *output)
{

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
	vtkMatrix4x4 *matrix;
	vtkPointData *pd, *outPD;

	TensorImageType::IndexType pixelIndex;
	TensorPixelType pixel;

	EigenValuesArrayType eigval;
	EigenVectorsMatrixType eigvec;

	double factor, norm_factor;
	float scalarValue;
	int i,j,k,p;
	int xSize, ySize, zSize;
	RealType cl, cp, cs;
	int numGlyphs = 0;
	double alpha, beta;

	vtkSuperquadricSource *superquadric;

	TensorImageType::PointType origin = this->input->GetOrigin();
	TensorImageType::SpacingType spacing = this->input->GetSpacing();


	// Generar la fuente de glifo escogida
	switch (this->GlyphType)
	{
		case ELLIPSOID: { 
			vtkSphereSource *sphere = vtkSphereSource::New();
			sphere->SetPhiResolution(this->PhiResolution);
			sphere->SetThetaResolution(this->ThetaResolution);
			this->source = sphere->GetOutput();
			this->source->Update();
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

	outPD = output->GetPointData();

	// Calcular el tamaño de la región a mostrar, limitada por los "Bounds"
	xSize = Bounds[1]-Bounds[0] + 1;
	ySize = Bounds[3]-Bounds[2] + 1;
	zSize = Bounds[5]-Bounds[4] + 1;
	numPts = xSize * ySize * zSize;

	sourcePts = source->GetPoints();
	numSourcePts = sourcePts->GetNumberOfPoints();
	numSourceCells = source->GetNumberOfCells();

	newPts = vtkPoints::New();

	pd = source->GetPointData();

	newScalars = vtkFloatArray::New();

	if ( (sourceNormals = pd->GetNormals()) )
	{
		newNormals = vtkFloatArray::New();
		newNormals->SetNumberOfComponents(3);
	}

	trans->PreMultiply();

	// Bucle principal, recorre la cuadrícula de la image
	// Se genera un glifo o ninguno en cada iteración
	for (i=this->Bounds[0]; i<=this->Bounds[1]; i++)
	for (j=this->Bounds[2]; j<=this->Bounds[3]; j++)
	for (k=this->Bounds[4]; k<=this->Bounds[5]; k++)
	{

		pixelIndex[0] = i;
		pixelIndex[1] = j;
		pixelIndex[2] = k;
		pixel = input->GetPixel(pixelIndex);

		pixel.ComputeEigenSystem(eigval,eigvec);

		// Si el tensor es nulo, no se representa
		if (eigval[0]==0) continue;

		pixel.ComputeShapeCoefficients(cl,cp,cs);

		// Discriminación de glifos
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

		// Comienza la transformación geométrica...

		trans->Identity();

		// Traslación del glifo a la posición adecuada
		x[0] = origin[0] + spacing[0] * pixelIndex[0];
		x[1] = origin[1] + spacing[1] * pixelIndex[1];
		x[2] = origin[2] + spacing[2] * pixelIndex[2];
		trans->Translate(x[0], x[1], x[2]);


		// Rotación del objeto en función de los autovectores
		// El factor signo evita la orientación incorrecta de las normales
		float signo = 1.0;
		if ( (eigvec[0][1]*eigvec[1][2]-eigvec[0][2]*eigvec[1][1]) * eigvec[2][0] < 0) signo = -1.0; 

		matrix->Element[0][0] = eigvec[0][0];
		matrix->Element[0][1] = eigvec[1][0];
		matrix->Element[0][2] = signo*eigvec[2][0];
		matrix->Element[1][0] = eigvec[0][1];
		matrix->Element[1][1] = eigvec[1][1];
		matrix->Element[1][2] = signo*eigvec[2][1];
		matrix->Element[2][0] = eigvec[0][2];
		matrix->Element[2][1] = eigvec[1][2];
		matrix->Element[2][2] = signo*eigvec[2][2];

		trans->Concatenate(matrix);

		// Cálculo del factor de normalización, el inverso del mayor autovalor
		if (fabs(eigval[0])>fabs(eigval[2]))
			norm_factor = 1 / eigval[0];

		else norm_factor = 1 / eigval[2];

		factor = norm_factor;

		if (this->Scaling)
			factor *= this->ScaleFactor;

		// Transformación del glifo en función de los autovalores
		trans->Scale(factor*eigval[0],factor*eigval[1],factor*eigval[2]);

 
		// Si el glifo es supercuádrico se modifican sus parámetros alfa y beta
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

		// Se aplica la transformación a los puntos y las normales
		trans->TransformPoints(sourcePts,newPts);

		if ( newNormals )
		{
			trans->TransformNormals(sourceNormals,newNormals);
		}

		// Cálculo del escalar para el coloreado
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

		for (p=0; p < numSourcePts; p++) 
		{
			newScalars->InsertNextValue(scalarValue);
		}

		// Se lleva la cuenta del número de glifos creado
		numGlyphs++;
	}


	if (this->inputPoints) numPts = this->inputPoints->GetNumberOfPoints();


	// Segundo bucle, para puntos de entrada específicos
	// El tensor se interpola, y el resto del proceso es idéntico
	if (this->inputPoints)

	for (inPtId=0; inPtId < numPts; inPtId++)
	{

		this->inputPoints->GetPoint(inPtId,x);
		interpolacionLogEuclidea(x,&pixel);

		pixel.ComputeEigenSystem(eigval,eigvec);

		if (eigval[0]==0) continue;

		pixel.ComputeShapeCoefficients(cl,cp,cs);

		trans->Identity();

		trans->Translate(x[0], x[1], x[2]);

		float signo = 1.0;
		if ( (eigvec[0][1]*eigvec[1][2]-eigvec[0][2]*eigvec[1][1]) * eigvec[2][0] < 0) signo = -1.0; 

		matrix->Element[0][0] = eigvec[0][0];
		matrix->Element[0][1] = eigvec[1][0];
		matrix->Element[0][2] = signo*eigvec[2][0];
		matrix->Element[1][0] = eigvec[0][1];
		matrix->Element[1][1] = eigvec[1][1];
		matrix->Element[1][2] = signo*eigvec[2][1];
		matrix->Element[2][0] = eigvec[0][2];
		matrix->Element[2][1] = eigvec[1][2];
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

		trans->TransformPoints(sourcePts,newPts);

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

		for (i=0; i < numSourcePts; i++) 
		{
			newScalars->InsertNextValue(scalarValue);
		}

		numGlyphs++;

	}

	// Se preparan las celdas del objeto de salida
	if ( (sourceCells=source->GetVerts())->GetNumberOfCells() > 0 )
	{
		cells = vtkCellArray::New();
		cells->Allocate(numGlyphs*sourceCells->GetSize());
		output->SetVerts(cells);
		cells->Delete();
	}
	if ( (sourceCells=source->GetLines())->GetNumberOfCells() > 0 )
	{
		cells = vtkCellArray::New();
		cells->Allocate(numGlyphs*sourceCells->GetSize());
		output->SetLines(cells);
		cells->Delete();
	}
	if ( (sourceCells=source->GetPolys())->GetNumberOfCells() > 0 )
	{
		cells = vtkCellArray::New();
		cells->Allocate(numGlyphs*sourceCells->GetSize());
		output->SetPolys(cells);
		cells->Delete();
	}
	if ( (sourceCells=source->GetStrips())->GetNumberOfCells() > 0 )
	{
		cells = vtkCellArray::New();
		cells->Allocate(numGlyphs*sourceCells->GetSize());
		output->SetStrips(cells);
		cells->Delete();
	}

	// Se copia la topología en el objeto de salida
	for (inPtId=0; inPtId < numGlyphs; inPtId++)
	{
		ptIncr = inPtId * numSourcePts;
		for (cellId=0; cellId < numSourceCells; cellId++)
		{
			cell = source->GetCell(cellId);
			cellPts = cell->GetPointIds();
			npts = cellPts->GetNumberOfIds();
			for (i=0; i < npts; i++)
			{
				pts[i] = cellPts->GetId(i) + ptIncr;
			}
			output->InsertNextCell(cell->GetCellType(),npts,pts);
		}
	}

	// Se actualizan la salida y se eliminan los objetos creados	

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

	return;
}

void vtkTensorGlyphDTI::interpolacionLineal (double x[3],TensorPixelType *pixel) {

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

void vtkTensorGlyphDTI::interpolacionLogEuclidea (double x[3],TensorPixelType *pixel) {

	TensorPixelType puntero;
	TensorPixelType pixel1, pixel2, pixel3, pixel4, pixel5, pixel6, pixel7, pixel8;
	TensorImageType::IndexType pixelIndex;
	TensorImageType::PointType origin = this->input->GetOrigin();
	TensorImageType::SpacingType spacing = this->input->GetSpacing();
	RealType t1,t2,t3;

	pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
	pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);
	pixelIndex[2] = (int) ((x[2]-origin[2]) / spacing[2]);

	t1 = ( (x[0]-origin[0]) / spacing[0] - pixelIndex[0] );
	t2 = ( (x[1]-origin[1]) / spacing[1] - pixelIndex[1] );
	t3 = ( (x[2]-origin[2]) / spacing[2] - pixelIndex[2] );

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

