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

//vtkCxxRevisionMacro(vtkSaturnTensorGlyph, "$Revision: 1.57.12.1 $");
//vtkStandardNewMacro(vtkSaturnTensorGlyph);

// Construct object with scaling on and scale factor 1.0. Eigenvalues are 
// extracted, glyphs are colored with input scalar data, and logarithmic
// scaling is turned off.
vtkSaturnTensorGlyph::vtkSaturnTensorGlyph()
{

  this->input = NULL;
  this->inputPoints = NULL;

  this->Scaling = 1;
  this->ScaleFactor = 1.0;
  this->ColorMode = COLOR_BY_SCALARS;
  this->ClampScaling = 0;
  this->MaxScaleFactor = 100;

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

  this->csThreshold = 1;
}

vtkSaturnTensorGlyph::~vtkSaturnTensorGlyph()
{
}

vtkPolyData *vtkSaturnTensorGlyph::GetOutput()
{

/*  // get the input and ouptut
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *source = vtkPolyData::SafeDownCast(
    sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
*/
//cout<<"1\n";
  vtkPolyData *output = vtkPolyData::New();
//cout<<"2\n";
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
  double factor;
  float scalarValue;
  int i,j,k;
  int xSize, ySize, zSize;
  RealType cl, cp, cs, t1, t2, t3, num;
  int cont=0;
  double alpha, beta;

  vtkSuperquadricSource *superquadric;

  TensorImageType::PointType origin = this->input->GetOrigin();
  TensorImageType::SpacingType spacing = this->input->GetSpacing();
//cout<<"1\n";
  switch (this->GlyphType)
    {
    case ELLIPSOID: { vtkSphereSource *sphere = vtkSphereSource::New();
      sphere->SetPhiResolution(this->PhiResolution);
      sphere->SetThetaResolution(this->ThetaResolution);
      this->source = sphere->GetOutput();
      this->source->Update();
//      sphere->Delete();
      break;}

    case CUBOID: {vtkCubeSource *cube = vtkCubeSource::New();
      source = cube->GetOutput();
      this->source->Update();
      break;}

    case SUPERQUADRIC: {superquadric = vtkSuperquadricSource::New();
      superquadric->SetPhiResolution(this->PhiResolution);
      superquadric->SetThetaResolution(this->ThetaResolution);
      this->source = superquadric->GetOutput();
      this->source->Update();
      break;}
    }

  //cout<<"puntos: "<<source->GetNumberOfCells()<<"\n";

  pts = new vtkIdType[source->GetMaxCellSize()];
  trans = vtkTransform::New();
  matrix = vtkMatrix4x4::New();
//cout<<"3\n";
  // set up working matrices
  m[0] = m0; m[1] = m1; m[2] = m2; 
  v[0] = v0; v[1] = v1; v[2] = v2; 

//  vtkDebugMacro(<<"Generating tensor glyphs");

  outPD = output->GetPointData();

  xSize = Bounds[1]-Bounds[0] + 1;
  ySize = Bounds[3]-Bounds[2] + 1;
  zSize = Bounds[5]-Bounds[4] + 1;
  numPts = xSize * ySize * zSize;
//cout<<"4\n";
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
//cout<<"7\n";
  if ( (sourceNormals = pd->GetNormals()) )
    {
    newNormals = vtkFloatArray::New();
    newNormals->SetNumberOfComponents(3);
    }


  //
  // Traverse all Input points, transforming glyph at Source points
  //
  trans->PreMultiply();
//cout<<"8\n";
  for (inPtId=0; inPtId < numPts; inPtId++)
    {
    ptIncr = inPtId * numSourcePts;

    // Translation is postponed

/*
    inTensors->GetTuple(inPtId, tensor);

    // compute orientation vectors and scale factors from tensor
    for (j=0; j<3; j++)
      {
      for (i=0; i<3; i++)
        {
        m[i][j] = tensor[i+3*j];
        }
      }
    vtkMath::Jacobi(m, w, v);

    //copy eigenvectors
    xv[0] = v[0][0]; xv[1] = v[1][0]; xv[2] = v[2][0];
    yv[0] = v[0][1]; yv[1] = v[1][1]; yv[2] = v[2][1];
    zv[0] = v[0][2]; zv[1] = v[1][2]; zv[2] = v[2][2];

    // compute scale factors
    w[0] *= this->ScaleFactor;
    w[1] *= this->ScaleFactor;
    w[2] *= this->ScaleFactor;
*/

    pixelIndex[0] = this->Bounds[0] + inPtId % xSize;
    pixelIndex[1] = this->Bounds[2] + inPtId / xSize % ySize;
    pixelIndex[2] = this->Bounds[4] + inPtId / (xSize * ySize) % zSize;

//cout<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<"\n";

    pixel = input->GetPixel(pixelIndex);

//    pixel.ComputeEigenSystem(eigval,eigvec);
/*
    m[0][0]=pixel[0];
    m[0][1]=pixel[1];
    m[0][2]=pixel[2];
    m[1][0]=pixel[1];
    m[1][1]=pixel[3];
    m[1][2]=pixel[4];
    m[2][0]=pixel[2];
    m[2][1]=pixel[4];
    m[2][2]=pixel[5];
    vtkMath::Jacobi(m, w, v);

cout<<v[0][0]<<"\t"<<v[1][0]<<"\t"<<v[2][0]<<"\t"<<v[0][1]<<"\t"<<v[1][1]<<"\t"<<v[2][1]<<"\n";
*/
pixel.ComputeEigenSystem(eigval,eigvec);
/*
cout<<pixel[0]<<"\t"<<pixel[1]<<"\t"<<pixel[2]<<"\t"<<pixel[3]<<"\t"<<pixel[4]<<"\t"<<pixel[5]<<"\n";
cout<<eigval[0]<<"\t"<<eigval[1]<<"\t"<<eigval[2]<<"\n";
cout<<eigvec[0][0]<<"\t"<<eigvec[1][0]<<"\t"<<eigvec[2][0]<<"\n";
cout<<eigvec[0][1]<<"\t"<<eigvec[1][1]<<"\t"<<eigvec[2][1]<<"\n";
cout<<eigvec[0][2]<<"\t"<<eigvec[1][2]<<"\t"<<eigvec[2][2]<<"\n\n";
*/

//cout<<w[0]<<" "<<w[1]<<" "<<w[2]<<"\t";
//cout<<v[0][0]<<" "<<v[0][1]<<" "<<v[0][2]<<"\t";
//cout<<v[1][0]<<" "<<v[1][1]<<" "<<v[1][2]<<"\t";
//cout<<v[2][0]<<" "<<v[2][1]<<" "<<v[2][2]<<"\t";
//cout<<atan2(v[1][0],v[0][0])*180/3.1415<<" "<<atan2(v[2][0],v[1][0])*180/3.1415<<" "<<atan2(v[0][0],v[2][0])*180/3.1415<<"\n";

    if (eigval[0]==0) continue;

    if (fabs(eigval[0])>fabs(eigval[2]))
      factor = 1 / eigval[0];

    else factor = 1 / eigval[2];
/*
    else if (w[0]==0)
      factor = 1;

    else factor = 1 / w[2];
*/
/*
    for (i=0; i<3; i++)
      {
        w[i] *= factor;
      }
*/



cout<<pixel[0]<<" "<<pixel[1]<<" "<<pixel[2]<<" "<<pixel[3]<<" "<<pixel[4]<<" "<<pixel[5]<<"\n";
cout<<eigval[0]<<" "<<eigval[1]<<" "<<eigval[2]<<"\t";
cout<<eigvec[0][0]<<" "<<eigvec[0][1]<<" "<<eigvec[0][2]<<"\t";
cout<<eigvec[1][0]<<" "<<eigvec[1][1]<<" "<<eigvec[1][2]<<"\t";
cout<<eigvec[2][0]<<" "<<eigvec[2][1]<<" "<<eigvec[2][2]<<"\n";




//    if ( (this->GlyphType==SUPERQUADRIC) || (this->ColorMode==COLOR_BY_CL) )
    pixel.ComputeShapeCoefficients(cl,cp,cs);

    if (fabs(cs)>=csThreshold) continue;
//if ( (w[0] < 0) || (w[1] < 0) || (w[2] < 0) )
//  cout<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<" "<<w[0]<<" "<<w[1]<<" "<<w[2]<<"\t"<<pixel[0]<<" "<<pixel[1]<<" "<<pixel[2]<<" "<<pixel[3]<<" "<<pixel[4]<<" "<<pixel[5]<<" "<<fabs(cs)<<"\n";


/* RECUPERAR ESTA PARTE SI SE SOLUCIONA EL CÁLCULO DE AUTOVECTORES
    if (eigval[0]!=0)
      factor = 1 / eigval[0];

    else factor=1;

    for (i=0; i<3; i++)
      {
        eigval[i] *= factor;
      }
*/

/*
    if ( this->ClampScaling )
      {
      for (maxScale=0.0, i=0; i<3; i++)
        {
        if ( maxScale < fabs(eigval[i]) )
          {
          maxScale = fabs(eigval[i]);
          }
        }
      if ( maxScale > this->MaxScaleFactor )
        {
        maxScale = this->MaxScaleFactor / maxScale;
        for (i=0; i<3; i++)
          {
          eigval[i] *= maxScale; //preserve overall shape of glyph
          }
        }
      }

    // normalization is postponed

    // make sure scale is okay (non-zero) and scale data
    for (maxScale=0.0, i=0; i<3; i++)
      {
      if ( eigval[i] > maxScale )
        {
        maxScale = eigval[i];
        }
      }
    if ( maxScale == 0.0 )
      {
      maxScale = 1.0;
      }
    for (i=0; i<3; i++)
      {
      if ( eigval[i] == 0.0 )
        {
        eigval[i] = maxScale * 1.0e-06;
        }
      }
*/
    // Now do the real work for each "direction"

    // Remove previous scales ...
    trans->Identity();

    // translate Source to Input point
//    input->GetPoint(inPtId, x);
    x[0] = origin[0] + spacing[0] * pixelIndex[0];
    x[1] = origin[1] + spacing[1] * pixelIndex[1];
    x[2] = origin[2] + spacing[2] * pixelIndex[2];
    trans->Translate(x[0], x[1], x[2]);

    // normalized eigenvectors rotate object for eigen direction 0

//trans->RotateX(atan2(1,1)*180/3.1415);
//trans->RotateZ(atan2(1,1)*180/3.1415);
//trans->RotateY(atan2(1,1)*180/3.1415);

//    trans->RotateX(atan2(v[1][0],v[0][0])*180/3.1415);
//    trans->RotateZ(atan2(v[1][0],v[0][0])*180/3.1415);
//    trans->RotateY(90+atan2(v[0][0],v[2][0])*180/3.1415);
//    trans->RotateX(90);
//    trans->RotateX(atan2(v[2][0],v[1][0])*180/3.1415);

    float signo = 1.0;
    if ( (eigvec[0][1]*eigvec[1][2]-eigvec[0][2]*eigvec[1][1]) * eigvec[2][0] < 0) signo = -1.0; 
/*
    matrix->Element[0][0] = v[0][0];
    matrix->Element[0][1] = v[0][1];
    matrix->Element[0][2] = v[0][2];
    matrix->Element[1][0] = v[1][0];
    matrix->Element[1][1] = v[1][1];
    matrix->Element[1][2] = v[1][2];
    matrix->Element[2][0] = v[2][0];
    matrix->Element[2][1] = v[2][1];
    matrix->Element[2][2] = v[2][2];
*/
/* RECUPERAR ESTA PARTE SI SE SOLUCIONA EL CÁLCULO DE AUTOVECTORES
*/
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

/* RECUPERAR ESTA PARTE SI SE SOLUCIONA EL CÁLCULO DE AUTOVECTORES
    trans->Scale(eigval[0], eigval[1], eigval[2]);
*/

    trans->Scale(signo*factor*eigval[0],factor*eigval[1],factor*eigval[2]);
//cout<<"a\n";
 
    if (this->GlyphType==SUPERQUADRIC) 
    {
      if (cl<=cp) {
        alpha = pow(1-cp,Gamma);
        beta = pow(1-cl,Gamma);
      }
      else {
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
//cout<<"gg\n";
    // Apply the transformation to a series of points, 
    // and append the results to outPts.
    if ( newNormals )
      {
      trans->TransformNormals(sourceNormals,newNormals);
      }

    switch (this->ColorMode)
      {

        case COLOR_BY_FA: scalarValue = pixel.GetFractionalAnisotropy();
          break;

        case COLOR_BY_RA: scalarValue = pixel.GetRelativeAnisotropy();
          break;

        case COLOR_BY_CL: scalarValue = cl;
          break;

      }

      // Copy point data from source
    for (i=0; i < numSourcePts; i++) 
      {
        newScalars->InsertNextValue(scalarValue);
//      newScalars->InsertTuple(ptIncr+i, &scalarValue);
//      newScalars->InsertTuple(cont*numSourcePts+i, &scalarValue);
      }

    cont++;
    }

  if (this->inputPoints) numPts = this->inputPoints->GetNumberOfPoints();

  if (this->inputPoints)

  for (inPtId=0; inPtId < numPts; inPtId++)
//  for (inPtId=0; inPtId < 1; inPtId++)
    {
    ptIncr = inPtId * numSourcePts;

/*
    pixelIndex[0] = this->Bounds[0] + inPtId % xSize;
    pixelIndex[1] = this->Bounds[2] + inPtId / xSize % ySize;
    pixelIndex[2] = this->Bounds[4] + inPtId / (xSize * ySize) % zSize;

    x[0] = origin[0] + spacing[0] * pixelIndex[0];
    x[1] = origin[1] + spacing[1] * pixelIndex[1];
    x[2] = origin[2] + spacing[2] * pixelIndex[2];
*/

    this->inputPoints->GetPoint(inPtId,x);
    
    pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
    pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);
    pixelIndex[2] = (int) ((x[2]-origin[2]) / spacing[2]);

pixel = input->GetPixel(pixelIndex);

    t1 = ( (x[0]-origin[0]) / spacing[0] - pixelIndex[0] ) / spacing[0];
    t2 = ( (x[1]-origin[1]) / spacing[1] - pixelIndex[1] ) / spacing[1];
    t3 = ( (x[2]-origin[2]) / spacing[2] - pixelIndex[2] ) / spacing[2];

/*
    pixel1 = (input->GetPixel(pixelIndex));
//    pixel1 = expTensor(pixel1);

    pixelIndex[0] = pixelIndex[0] + 1;

    pixel2 = (input->GetPixel(pixelIndex));
//    pixel2 = expTensor(pixel2);

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] + 1;

    pixel3 = (input->GetPixel(pixelIndex));
//    pixel3 = expTensor(pixel3);

    pixelIndex[0] = pixelIndex[0] + 1;

    pixel4 = (input->GetPixel(pixelIndex));
//    pixel4 = expTensor(pixel4);

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] - 1;
    pixelIndex[2] = pixelIndex[2] + 1;

    pixel5 = (input->GetPixel(pixelIndex));
//    pixel5 = expTensor(pixel5);

    pixelIndex[0] = pixelIndex[0] + 1;

    pixel6 = (input->GetPixel(pixelIndex));
//    pixel6 = expTensor(pixel6);

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] + 1;

    pixel7 = (input->GetPixel(pixelIndex));
//    pixel7 = expTensor(pixel7);

    pixelIndex[0] = pixelIndex[0] + 1;

    pixel8 = (input->GetPixel(pixelIndex));
//    pixel8 = expTensor(pixel8);

//    pixel = (pixel1*realValue + pixel2*(1-realValue)) * realValue + (pixel3*realValue + pixel4*(1-realValue)) * (1-realValue);
    pixel = ((pixel1*t1 + pixel2*(1-t1)) * t2 + (pixel3*t1 + pixel4*(1-t1)) * (1-t2)) * t3 + ((pixel5*t1 + pixel6*(1-t1)) * t2 + (pixel7*t1 + pixel8*(1-t1)) * (1-t2)) * (1-t3);
//    pixel = logTensor(pixel);
//cout<<"todos los pixels 548 "<<pixel1<<"\n"<<pixel2<<"\n"<<pixel3<<"\n"<<pixel4<<"\n"<<pixel5<<"\n"<<pixel6<<"\n"<<pixel7<<"\n"<<pixel8<<"\n\n";
*/
    interpolacionLogEuclidea(x,&pixel);
//interpolacionLineal(x,&pixel);
//cout<<pixel<<"\n";
 
//cout<<"e\n";
/*
    m[0][0]=pixel[0];
    m[0][1]=pixel[1];
    m[0][2]=pixel[2];
    m[1][0]=pixel[1];
    m[1][1]=pixel[3];
    m[1][2]=pixel[4];
    m[2][0]=pixel[2];
    m[2][1]=pixel[4];
    m[2][2]=pixel[5];
    vtkMath::Jacobi(m, w, v);
*/
pixel.ComputeEigenSystem(eigval,eigvec);

/*
if (fabs(v[0][0])-fabs(eigvec[0][0]) > 0.0001) {
cout<<w[0]<<"\t"<<eigval[0]<<"\t"<<w[1]<<"\t"<<eigval[1]<<"\t"<<w[2]<<"\t"<<eigval[2]<<"\n";
cout<<v[0][0]<<"\t"<<eigvec[0][0]<<"\t"<<v[1][0]<<"\t"<<eigvec[1][0]<<"\t"<<v[2][0]<<"\t"<<eigvec[2][0]<<"\n";
cout<<v[0][1]<<"\t"<<eigvec[0][1]<<"\t"<<v[1][1]<<"\t"<<eigvec[1][1]<<"\t"<<v[2][1]<<"\t"<<eigvec[2][1]<<"\n";
cout<<v[0][2]<<"\t"<<eigvec[0][2]<<"\t"<<v[1][2]<<"\t"<<eigvec[1][2]<<"\t"<<v[2][2]<<"\t"<<eigvec[2][2]<<"\n\n";
}
*/
    if (eigval[0]==0) continue;

    if (fabs(eigval[0])>fabs(eigval[2]))
      factor = 1 / eigval[0];

    else factor = 1 / eigval[2];

    pixel.ComputeShapeCoefficients(cl,cp,cs);

    if (fabs(cs)>=csThreshold) continue;

    // Now do the real work for each "direction"

    // Remove previous scales ...
    trans->Identity();

    // translate Source to Input point
/*
    x[0] = origin[0] + spacing[0] * pixelIndex[0];
    x[1] = origin[1] + spacing[1] * pixelIndex[1] - 0.5 * spacing[1];
    x[2] = origin[2] + spacing[2] * pixelIndex[2] - 0.5 * spacing[2];
*/
    trans->Translate(x[0], x[1], x[2]);

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

    float signo = 1.0;
    if ( (eigvec[0][1]*eigvec[1][2]-eigvec[0][2]*eigvec[1][1]) * eigvec[2][0] < 0) signo = -1.0; 

    trans->Scale(signo*factor*eigval[0],factor*eigval[1],factor*eigval[2]);
//cout<<"ee\n";
    if (this->GlyphType==SUPERQUADRIC) 
    {
      if (cl<=cp) {
        alpha = pow(1-cp,Gamma);
        beta = pow(1-cl,Gamma);
      }
      else {
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
//cout<<"ff\n";
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

        case COLOR_BY_FA: scalarValue = pixel.GetFractionalAnisotropy();
          break;

        case COLOR_BY_RA: scalarValue = pixel.GetRelativeAnisotropy();
          break;

        case COLOR_BY_CL: scalarValue = cl;
          break;

      }

      // Copy point data from source
    for (i=0; i < numSourcePts; i++) 
      {
        newScalars->InsertNextValue(scalarValue);
      }
    cont++;
    }

//  vtkDebugMacro(<<"Generated " << numPts <<" tensor glyphs");
//cout<<"gg\n";
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
//cout<<"hh\n";
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
//cout<<"ii\n";
  //
  // Update output and release memory
  //
  delete [] pts;

  output->SetPoints(newPts);
  newPts->Delete();

  if ( newScalars )
    {
//    int idx = outPD->AddArray(newScalars);
//    outPD->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    output->GetPointData()->SetScalars(newScalars);
//    outPD->SetScalars(newScalars);
    newScalars->Delete();
    }

  if ( newNormals )
    {
    outPD->SetNormals(newNormals);
    newNormals->Delete();
    }


  trans->Delete();
  matrix->Delete();
//  source->Delete();
//cout<<"jj\n";
  return output;
}

void vtkSaturnTensorGlyph::expTensor(TensorPixelType *inPixel) {

//  TensorImageType::IndexType pixelIndex;
//  TensorPixelType pixel = input->GetPixel(pixelIndex);
 
//  EigenValuesArrayType eigval;
//  EigenVectorsMatrixType eigvec;
//  ComputeEigenSystem(eigval, eigvec);
  
  TensorPixelType pixel = *inPixel;
  
  double *m[3], w[3], *v[3];
  double m0[3], m1[3], m2[3];
  double v0[3], v1[3], v2[3];
  double xv[3], yv[3], zv[3];

  m[0] = m0; m[1] = m1; m[2] = m2; 
  v[0] = v0; v[1] = v1; v[2] = v2; 

  m[0][0]=pixel[0];
  m[0][1]=pixel[1];
  m[0][2]=pixel[2];
  m[1][0]=pixel[1];
  m[1][1]=pixel[3];
  m[1][2]=pixel[4];
  m[2][0]=pixel[2];
  m[2][1]=pixel[4];
  m[2][2]=pixel[5];
  vtkMath::Jacobi(m, w, v);

  double eigval[3];
	
  eigval[0]=exp(w[0]);
  eigval[1]=exp(w[1]);
  eigval[2]=exp(w[2]);

//  itk::Matrix<RealType,3,3> tensorMatrix;
//  tensorMatrix=logEigenVectors.GetVnlMatrix()*eigval.GetVnlMatrix()*logEigenVectors.GetTranspose();
//  tensorMatrix=logEigenVectors.GetTranspose()*eigval.GetVnlMatrix()*logEigenVectors.GetVnlMatrix();
//  tensorMatrix=logEigenVectors.GetTranspose()*eigval*logEigenVectors;
//v[0][0]*eigval[0];

  TensorPixelType result;	
  result[0]=v[0][0]*v[0][0]*eigval[0]+v[0][1]*v[0][1]*eigval[1]+v[0][2]*v[0][2]*eigval[2];
  result[1]=v[0][0]*v[1][0]*eigval[0]+v[0][1]*v[1][1]*eigval[1]+v[0][2]*v[1][2]*eigval[2];
  result[2]=v[0][0]*v[2][0]*eigval[0]+v[0][1]*v[2][1]*eigval[1]+v[0][2]*v[2][2]*eigval[2];
  result[3]=v[1][0]*v[1][0]*eigval[0]+v[1][1]*v[1][1]*eigval[1]+v[1][2]*v[1][2]*eigval[2];
  result[4]=v[1][0]*v[2][0]*eigval[0]+v[1][1]*v[2][1]*eigval[1]+v[1][2]*v[2][2]*eigval[2];
  result[5]=v[2][0]*v[2][0]*eigval[0]+v[2][1]*v[2][1]*eigval[1]+v[2][2]*v[2][2]*eigval[2];

  *inPixel = result;

}

void vtkSaturnTensorGlyph::logTensor(TensorPixelType *inPixel) {

//  TensorImageType::IndexType pixelIndex;
//  TensorPixelType pixel = input->GetPixel(pixelIndex);
 
//  EigenValuesArrayType eigval;
//  EigenVectorsMatrixType eigvec;
//  ComputeEigenSystem(eigval, eigvec);

  TensorPixelType pixel = *inPixel;

  double *m[3], w[3], *v[3];
  double m0[3], m1[3], m2[3];
  double v0[3], v1[3], v2[3];
  double xv[3], yv[3], zv[3];

  m[0] = m0; m[1] = m1; m[2] = m2; 
  v[0] = v0; v[1] = v1; v[2] = v2; 

  m[0][0]=pixel[0];
  m[0][1]=pixel[1];
  m[0][2]=pixel[2];
  m[1][0]=pixel[1];
  m[1][1]=pixel[3];
  m[1][2]=pixel[4];
  m[2][0]=pixel[2];
  m[2][1]=pixel[4];
  m[2][2]=pixel[5];
  vtkMath::Jacobi(m, w, v);

  double eigval[3];
	
  eigval[0]=log(w[0]);
  eigval[1]=log(w[1]);
  eigval[2]=log(w[2]);

//  itk::Matrix<RealType,3,3> tensorMatrix;
//  tensorMatrix=logEigenVectors.GetVnlMatrix()*eigval.GetVnlMatrix()*logEigenVectors.GetTranspose();
//  tensorMatrix=logEigenVectors.GetTranspose()*eigval.GetVnlMatrix()*logEigenVectors.GetVnlMatrix();
//  tensorMatrix=logEigenVectors.GetTranspose()*eigval*logEigenVectors;
//v[0][0]*eigval[0];

  TensorPixelType result;	
  result[0]=v[0][0]*v[0][0]*eigval[0]+v[0][1]*v[0][1]*eigval[1]+v[0][2]*v[0][2]*eigval[2];
  result[1]=v[0][0]*v[1][0]*eigval[0]+v[0][1]*v[1][1]*eigval[1]+v[0][2]*v[1][2]*eigval[2];
  result[2]=v[0][0]*v[2][0]*eigval[0]+v[0][1]*v[2][1]*eigval[1]+v[0][2]*v[2][2]*eigval[2];
  result[3]=v[1][0]*v[1][0]*eigval[0]+v[1][1]*v[1][1]*eigval[1]+v[1][2]*v[1][2]*eigval[2];
  result[4]=v[1][0]*v[2][0]*eigval[0]+v[1][1]*v[2][1]*eigval[1]+v[1][2]*v[2][2]*eigval[2];
  result[5]=v[2][0]*v[2][0]*eigval[0]+v[2][1]*v[2][1]*eigval[1]+v[2][2]*v[2][2]*eigval[2];
	
  *inPixel = result;

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

EigenValuesArrayType eigval;
EigenVectorsMatrixType eigvec;
    
    pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
    pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);
    pixelIndex[2] = (int) ((x[2]-origin[2]) / spacing[2]);

    t1 = ( (x[0]-origin[0]) / spacing[0] - pixelIndex[0] ) / spacing[0];
    t2 = ( (x[1]-origin[1]) / spacing[1] - pixelIndex[1] ) / spacing[1];
    t3 = ( (x[2]-origin[2]) / spacing[2] - pixelIndex[2] ) / spacing[2];

    pixel1 = (input->GetPixel(pixelIndex));
pixel1.ComputeEigenSystem(eigval,eigvec);
cout<<"Empezamos:\n"<<pixel1<<"\n";
cout<<eigval<<"\n"<<eigvec<<"\n";
    pixel1 = pixel1.ComputeTensorFromLogVector();
pixel1.ComputeEigenSystem(eigval,eigvec);
cout<<pixel1<<"\n";
cout<<eigval<<"\n"<<eigvec<<"\n";

    pixelIndex[0] = pixelIndex[0] + 1;
    pixel2 = (input->GetPixel(pixelIndex));
    pixel2 = pixel2.ComputeTensorFromLogVector();
cout<<pixel2<<"\n";

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] + 1;
    pixel3 = (input->GetPixel(pixelIndex));
    pixel3 = pixel3.ComputeTensorFromLogVector();
cout<<pixel3<<"\n";

    pixelIndex[0] = pixelIndex[0] + 1;
    pixel4 = (input->GetPixel(pixelIndex));
    pixel4 = pixel4.ComputeTensorFromLogVector();
cout<<pixel4<<"\n";

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] - 1;
    pixelIndex[2] = pixelIndex[2] + 1;
    pixel5 = (input->GetPixel(pixelIndex));
    pixel5 = pixel5.ComputeTensorFromLogVector();
cout<<pixel5<<"\n";

    pixelIndex[0] = pixelIndex[0] + 1;
    pixel6 = (input->GetPixel(pixelIndex));
    pixel6 = pixel6.ComputeTensorFromLogVector();
cout<<pixel6<<"\n";

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] + 1;
    pixel7 = (input->GetPixel(pixelIndex));
    pixel7 = pixel7.ComputeTensorFromLogVector();
cout<<pixel7<<"\n";

    pixelIndex[0] = pixelIndex[0] + 1;
    pixel8 = (input->GetPixel(pixelIndex));
    pixel8 = pixel8.ComputeTensorFromLogVector();
//cout<<pixel8<<"\n"<<&puntero<<"\n";
//    expTensor(&pixel8);
//cout<<&pixel8<<"\n";
cout<<pixel1<<"\n";
cout<<pixel2<<"\n";
cout<<pixel3<<"\n";
cout<<pixel4<<"\n";
cout<<pixel5<<"\n";
cout<<pixel6<<"\n";
cout<<pixel7<<"\n";
cout<<pixel8<<"\n";

    *pixel = ((pixel1*(1-t1) + pixel2*t1) * (1-t2) + (pixel3*(1-t1) + pixel4*t1) * t2) * (1-t3) + ((pixel5*(1-t1) + pixel6*t1) * (1-t2) + (pixel7*(1-t1) + pixel8*t1) * t2) * t3;
pixel->ComputeEigenSystem(eigval,eigvec);
cout<<*pixel<<"\n";
cout<<eigval<<"\n"<<eigvec<<"\n";

*pixel = pixel->ComputeLogVector();
pixel->ComputeEigenSystem(eigval,eigvec);
cout<<*pixel<<"\n";
cout<<eigval<<"\n"<<eigvec<<"\n";

//    logTensor(pixel);
//cout<<*pixel<<"\n\n";
 
}

/*
void vtkSaturnTensorGlyph::interpolacionLogEuclidea (double x[3],TensorPixelType *pixel_in) {

    TensorPixelType *puntero, *pixel;
    TensorPixelType *pixel1, *pixel2, *pixel3, *pixel4, *pixel5, *pixel6, *pixel7, *pixel8;
    TensorImageType::IndexType pixelIndex;
    TensorImageType::PointType origin = this->input->GetOrigin();
    TensorImageType::SpacingType spacing = this->input->GetSpacing();
    RealType t1,t2,t3;

EigenValuesArrayType eigval;
EigenVectorsMatrixType eigvec;
    
    pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
    pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);
    pixelIndex[2] = (int) ((x[2]-origin[2]) / spacing[2]);

    t1 = ( (x[0]-origin[0]) / spacing[0] - pixelIndex[0] ) / spacing[0];
    t2 = ( (x[1]-origin[1]) / spacing[1] - pixelIndex[1] ) / spacing[1];
    t3 = ( (x[2]-origin[2]) / spacing[2] - pixelIndex[2] ) / spacing[2];

cout<<"t1,t2,t3 "<<t1<<"\t"<<t2<<"\t"<<t3<<"\n";

    pixel = &(input->GetPixel(pixelIndex));
    pixel->ComputeEigenSystem(eigval,eigvec);
//cout<<"Aqui llega ";
//cout<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n";
cout<<"pixel 1 a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

    pixel1 = new TensorPixelType(pixel->ComputeTensorFromLogVector());
    pixel1->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 1 b "<<*pixel1<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

//pixel1->ComputeEigenSystem(eigval,eigvec);
//cout<<pixel1[0]<<" "<<pixel1[1]<<" "<<pixel1[2]<<" "<<pixel1[3]<<" "<<pixel1[4]<<" "<<pixel1[5]<<"\n";
//cout<<eigval[0]<<" "<<eigval[1]<<" "<<eigval[2]<<"\t";
//cout<<eigvec[0][0]<<" "<<eigvec[0][1]<<" "<<eigvec[0][2]<<"\t";
//cout<<eigvec[1][0]<<" "<<eigvec[1][1]<<" "<<eigvec[1][2]<<"\t";
//cout<<eigvec[2][0]<<" "<<eigvec[2][1]<<" "<<eigvec[2][2]<<"\n";
//cout<<"hola?? "<<pixel1[0]<<" "<<pixel1[1]<<" "<<pixel1[2]<<" "<<pixel1[3]<<" "<<pixel1[4]<<" "<<pixel1[5]<<"\n";

    pixelIndex[0] = pixelIndex[0] + 1;

    pixel = &(input->GetPixel(pixelIndex));
    pixel->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 2 a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";
    pixel2 = new TensorPixelType(pixel->ComputeTensorFromLogVector());
    pixel2->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 2 b "<<*pixel2<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] + 1;

    pixel = &(input->GetPixel(pixelIndex));
    pixel->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 3 a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";
    pixel3 = new TensorPixelType(pixel->ComputeTensorFromLogVector());
    pixel3->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 3 b "<<*pixel3<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

    pixelIndex[0] = pixelIndex[0] + 1;

    pixel = &(input->GetPixel(pixelIndex));
    pixel->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 4 a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";
    pixel4 = new TensorPixelType(pixel->ComputeTensorFromLogVector());
    pixel4->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 4 b "<<*pixel4<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] - 1;
    pixelIndex[2] = pixelIndex[2] + 1;

    pixel = &(input->GetPixel(pixelIndex));
    pixel->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 5 a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";
    pixel5 = new TensorPixelType(pixel->ComputeTensorFromLogVector());
    pixel5->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 5 b "<<*pixel5<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

    pixelIndex[0] = pixelIndex[0] + 1;

    pixel = &(input->GetPixel(pixelIndex));
    pixel->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 6 a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";
    pixel6 = new TensorPixelType(pixel->ComputeTensorFromLogVector());
    pixel6->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 6 b "<<*pixel6<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

    pixelIndex[0] = pixelIndex[0] - 1;
    pixelIndex[1] = pixelIndex[1] + 1;

    pixel = &(input->GetPixel(pixelIndex));
    pixel->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 7 a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";
    pixel7 = new TensorPixelType(pixel->ComputeTensorFromLogVector());
    pixel7->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 7 b "<<*pixel7<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

    pixelIndex[0] = pixelIndex[0] + 1;

    pixel = &(input->GetPixel(pixelIndex));
    pixel->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 8 a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";
    pixel8 = new TensorPixelType(pixel->ComputeTensorFromLogVector());
    pixel8->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel 8 b "<<*pixel8<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

cout<<"todos los pixels "<<*pixel1<<"\n"<<*pixel2<<"\n"<<*pixel3<<"\n"<<*pixel4<<"\n"<<*pixel5<<"\n"<<*pixel6<<"\n"<<*pixel7<<"\n"<<*pixel8<<"\n\n";


    *pixel = (((*pixel1)*(1-t1) + (*pixel2)*t1) * (1-t2) + ((*pixel3)*(1-t1) + (*pixel4)*t1) * t2) * (1-t3) + (((*pixel5)*(1-t1) + (*pixel6)*t1) * (1-t2) + ((*pixel7)*(1-t1) + (*pixel8)*t1) * t2) * t3;
pixel->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel out a "<<*pixel<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";

    *pixel_in = pixel->ComputeLogVector();
//cout<<*pixel<<"\n"<<*pixel_in<<"\n\n";
pixel_in->ComputeEigenSystem(eigval,eigvec);
cout<<"pixel out b "<<*pixel_in<<"\n"<<eigval<<"\t"<<eigvec<<"\n\n";
    
}
*/










