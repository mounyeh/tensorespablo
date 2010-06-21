/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkStrainTensorGlyph.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkStrainTensorGlyph.h"

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

#include <algorithm>

//vtkCxxRevisionMacro(vtkStrainTensorGlyph, "$Revision: 1.57.12.1 $");
//vtkStandardNewMacro(vtkStrainTensorGlyph);

// Construct object with scaling on and scale factor 1.0. Eigenvalues are 
// extracted, glyphs are colored with input scalar data, and logarithmic
// scaling is turned off.
vtkStrainTensorGlyph::vtkStrainTensorGlyph()
{
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

vtkStrainTensorGlyph::~vtkStrainTensorGlyph()
{
}

vtkPolyData *vtkStrainTensorGlyph::GetOutput()
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
  vtkPolyData *source = this->source;
  vtkPolyData *output = vtkPolyData::New();
  source->Update();
  input->Update();
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
//  double *m[3], w[3], *v[3];
//  double m0[3], m1[3], m2[3];
//  double v0[3], v1[3], v2[3];
//  double xv[3], yv[3], zv[3];
  double *m[2], w[2], *v[2];
  double m0[2], m1[2];
  double v0[2], v1[2];
  double xv[2], yv[2], zv[2];
  double maxScale;
  vtkPointData *pd, *outPD;

//  TensorImageType::IndexType pixelIndex;
//  TensorPixelType pixel;
  EigenValuesArrayType eigval;
  EigenVectorsMatrixType eigvec;
  double factor;
  float scalarValue;
  float scalar2;
  int i,j,k;
  int xSize, ySize, zSize;
  RealType cl, cp, cs;
  int cont=0;
  double def;
  double angle;

//  TensorImageType::PointType origin = this->input->GetOrigin();
//  TensorImageType::SpacingType spacing = this->input->GetSpacing();

//  numDirs = (this->ThreeGlyphs?3:1)*(this->Symmetric+1);

  pts = new vtkIdType[source->GetMaxCellSize()];
  trans = vtkTransform::New();
  matrix = vtkMatrix4x4::New();
//cout<<"3\n";
  // set up working matrices
  m[0] = m0; m[1] = m1; //m[2] = m2; 
  v[0] = v0; v[1] = v1; //v[2] = v2; 

//  vtkDebugMacro(<<"Generating tensor glyphs");

  pd = input->GetPointData();
  outPD = output->GetPointData();
  inTensors = pd->GetTensors();
  inScalars = pd->GetScalars();
  numPts = input->GetNumberOfPoints();

//  xSize = Bounds[1]-Bounds[0] + 1;
//  ySize = Bounds[3]-Bounds[2] + 1;
//  zSize = Bounds[5]-Bounds[4] + 1;
//  numPts = xSize * ySize * zSize;
//  numPts = xSize * ySize * zSize + 1;
//cout<<"4\n";
  //
  // Allocate storage for output PolyData
  //
  sourcePts = source->GetPoints();
  numSourcePts = sourcePts->GetNumberOfPoints();
  numSourceCells = source->GetNumberOfCells();

  newPts = vtkPoints::New();
//  newPts->Allocate(numPts*numSourcePts);
//  for (i=0; i<numSourcePts; i++) {
//    newPts->InsertNextPoint(sourcePts->GetPoint(i));
//  }

  // only copy scalar data through
  pd = source->GetPointData();

  newScalars = vtkFloatArray::New();
//  newScalars->Allocate(numPts*numSourcePts);
//cout<<"7\n";
  if ( (sourceNormals = pd->GetNormals()) )
    {
    newNormals = vtkFloatArray::New();
    newNormals->SetNumberOfComponents(3);
//    newNormals->Allocate(3*numPts*numSourcePts);
    }
  if ( newScalars )
    {
    int idx = outPD->AddArray(newScalars);
    outPD->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    newScalars->Delete();
    }

  if ( newNormals )
    {
    outPD->SetNormals(newNormals);
    newNormals->Delete();
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

    inTensors->GetTuple(inPtId, tensor);
//cout<<"a\n";

    scalarValue = inScalars->GetTuple1(inPtId);
/*
    // compute orientation vectors and scale factors from tensor
    for (j=0; j<2; j++)
      {
      for (i=0; i<2; i++)
        {
        m[i][j] = tensor[i+3*j];
        }
      }
//cout<<"b\n";
*/
//    def = m[0][0];
//cout<<"c\n";
/*
    m[0][0] = tensor[0];
    m[0][1] = tensor[1];
    m[1][0] = tensor[2];
    m[1][1] = tensor[3];

    vtkMath::JacobiN(m, 2, w, v);
*/

ComputeEigenSystem(inPtId,v,w);

//if (tensor[0]!=0)
//cout<<v[0][0]<<" "<<v[1][0]<<"\t\t"<<v[0][1]<<" "<<v[1][1]<<"\t\t"<<w[0]<<" "<<w[1]<<" "<<"\n";


    //copy eigenvectors
    xv[0] = v[0][0]; xv[1] = v[0][1];// xv[2] = v[2][0];
    yv[0] = v[1][0]; yv[1] = v[1][1];// yv[2] = v[2][1];
//    zv[0] = v[0][2]; zv[1] = v[1][2]; zv[2] = v[2][2];

//if (tensor[0]!=0)
//cout<<xv[0]<<" "<<yv[0]<<"\t\t"<<xv[1]<<" "<<yv[1]<<"\t\t"<<w[0]<<" "<<w[1]<<" "<<"\n";

    // compute scale factors
    w[0] *= this->ScaleFactor;
    w[1] *= this->ScaleFactor;
//    w[2] *= this->ScaleFactor;


//    pixelIndex[0] = this->Bounds[0] + inPtId % xSize;
//    pixelIndex[1] = this->Bounds[2] + inPtId / xSize % ySize;
//    pixelIndex[2] = this->Bounds[4] + inPtId / (xSize * ySize) % zSize;

//cout<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<"\n";

//    pixel = input->GetPixel(pixelIndex);

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
*/

    if (scalarValue==0) continue;

/*
    if (w[0]==0) continue;

*/
    if (fabs(w[0])>fabs(w[2]))
      factor = 1 / w[0];

    else factor = 1 / w[2];

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


/*
cout<<pixel[0]<<" "<<pixel[1]<<" "<<pixel[2]<<" "<<pixel[3]<<" "<<pixel[4]<<" "<<pixel[5]<<"\n";
cout<<eigval[0]<<" "<<eigval[1]<<" "<<eigval[2]<<"\t";
cout<<eigvec[0][0]<<" "<<eigvec[0][1]<<" "<<eigvec[0][2]<<"\t";
cout<<eigvec[1][0]<<" "<<eigvec[1][1]<<" "<<eigvec[1][2]<<"\t";
cout<<eigvec[2][0]<<" "<<eigvec[2][1]<<" "<<eigvec[2][2]<<"\n";
cout<<w[0]<<" "<<w[1]<<" "<<w[2]<<"\t";
cout<<v[0][0]<<" "<<v[0][1]<<" "<<v[0][2]<<"\t";
cout<<v[1][0]<<" "<<v[1][1]<<" "<<v[1][2]<<"\t";
cout<<v[2][0]<<" "<<v[2][1]<<" "<<v[2][2]<<"\n\n";
*/

//    if ( (this->GlyphType==SUPERQUADRIC) || (this->ColorMode==COLOR_BY_CL) )
/*
    pixel.ComputeShapeCoefficients(cl,cp,cs);

    if (fabs(cs)>=csThreshold) continue;
*/
//if ( (w[0] < 0) || (w[1] < 0) || (w[2] < 0) )
//  cout<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<" "<<w[0]<<" "<<w[1]<<" "<<w[2]<<"\t"<<pixel[0]<<" "<<pixel[1]<<" "<<pixel[2]<<" "<<pixel[3]<<" "<<pixel[4]<<" "<<pixel[5]<<"\n";


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
    input->GetPoint(inPtId, x);
/*
    x[0] = origin[0] + spacing[0] * pixelIndex[0];
    x[1] = origin[1] + spacing[1] * pixelIndex[1];
    x[2] = origin[2] + spacing[2] * pixelIndex[2];
*/
    trans->Translate(x[0], x[1], x[2]);

    // normalized eigenvectors rotate object for eigen direction 0
 
    double ang1 = atan2(yv[1],xv[1]) * 180 / 3.1415;
    double ang2 = atan(-yv[1]/fabs(xv[1])) * 180 / 3.1415;

//cout<<"Angulo: "<<atan2(yv[0],xv[0])*180/3.1415<<" "<<atan2(yv[0],-xv[0])*180/3.1415<<"\n";
//cout<<"Angulo: "<<ang1<<"\n";

//    trans->RotateX(ang1);
    trans->RotateY(ang1);

    double signo = -1.0;
    if (v[0][1]<0 && v[1][1]<0) signo = 1.0;
/*
    matrix->Element[0][0] = xv[1];
    matrix->Element[0][1] = 0;		//yv[0];
    matrix->Element[0][2] = yv[1];	//zv[0];
    matrix->Element[1][0] = 0;
    matrix->Element[1][1] = 1;
    matrix->Element[1][2] = 0;
    matrix->Element[2][0] = signo*xv[0];
    matrix->Element[2][1] = 0;		//yv[1];
    matrix->Element[2][2] = signo*yv[0];//zv[1];

    matrix->Element[0][0] = -1;
    matrix->Element[0][1] = 0;		//yv[0];
    matrix->Element[0][2] = -1;	//zv[0];
    matrix->Element[1][0] = 0;
    matrix->Element[1][1] = 1;
    matrix->Element[1][2] = 0;
    matrix->Element[2][0] = 1;
    matrix->Element[2][1] = 0;		//yv[1];
    matrix->Element[2][2] = -1;//zv[1];
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Ver página 64 del pdf
*/
/*
    angle = atan2(-x[0],-x[2]);
    trans->RotateY(90 + angle*180/3.1415);
*/

/* RECUPERAR ESTA PARTE SI SE SOLUCIONA EL CÁLCULO DE AUTOVECTORES
    matrix->Element[0][0] = eigvec[0][0];
    matrix->Element[0][1] = eigvec[1][0];
    matrix->Element[0][2] = eigvec[2][0];
    matrix->Element[1][0] = eigvec[0][1];
    matrix->Element[1][1] = eigvec[1][1];
    matrix->Element[1][2] = eigvec[2][1];
    matrix->Element[2][0] = eigvec[0][2];
    matrix->Element[2][1] = eigvec[1][2];
    matrix->Element[2][2] = eigvec[2][2];
*/
//    trans->Concatenate(matrix);

/* RECUPERAR ESTA PARTE SI SE SOLUCIONA EL CÁLCULO DE AUTOVECTORES
    trans->Scale(eigval[0], eigval[1], eigval[2]);
*/

//    trans->RotateY(90.0);

//    trans->Scale(factor*fabs(w[0]),factor*fabs(w[2]),.001);

//      trans->Scale(w[0],1,fabs(w[1]));
if (fabs(w[0])>fabs(w[1])) factor = w[0];
  else factor = w[1];
 
      trans->Scale(w[0]/factor,.2,fabs(w[1]/factor));

/*
    if (GlyphType==SUPERQUADRIC) 
      {
      source = this->GetSource();
      sourcePts = source->GetPoints();
      }
*/
    // multiply points (and normals if available) by resulting
    // matrix
    trans->TransformPoints(sourcePts,newPts);
    // Apply the transformation to a series of points, 
    // and append the results to outPts.
    if ( newNormals )
      {
      trans->TransformNormals(sourceNormals,newNormals);
      }

/*
    switch (this->ColorMode)
      {

        case COLOR_BY_FA: scalarValue = pixel.GetFractionalAnisotropy();
          break;

        case COLOR_BY_RA: scalarValue = pixel.GetRelativeAnisotropy();
          break;

        case COLOR_BY_CL: scalarValue = cl;
          break;

      }
*/

if (fabs(w[0])>fabs(w[1])) scalarValue = w[0];
  else scalarValue = w[1];


      // Copy point data from source
    for (i=0; i < numSourcePts; i++) 
      {
        newScalars->InsertNextTupleValue(&scalarValue);
//      newScalars->InsertTuple(ptIncr+i, &scalarValue);
//      newScalars->InsertTuple(cont*numSourcePts+i, &scalarValue);
      }
    cont++;
    }
//  vtkDebugMacro(<<"Generated " << numPts <<" tensor glyphs");

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

  trans->Delete();
  matrix->Delete();
//  source->Delete();

  return output;
}

vtkPolyData *vtkStrainTensorGlyph::GetDeformOutput()
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
  vtkPolyData *source = this->source;
  vtkPolyData *output = vtkPolyData::New();
  source->Update();
  input->Update();
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
//  double *m[3], w[3], *v[3];
//  double m0[3], m1[3], m2[3];
//  double v0[3], v1[3], v2[3];
//  double xv[3], yv[3], zv[3];
  double *m[2], w[2], *v[2];
  double m0[2], m1[2];
  double v0[2], v1[2];
  double xv[2], yv[2], zv[2];
  double maxScale;
  vtkPointData *pd, *outPD;

//  TensorImageType::IndexType pixelIndex;
//  TensorPixelType pixel;
  EigenValuesArrayType eigval;
  EigenVectorsMatrixType eigvec;
  double factor;
  float scalarValue;
  float scalar2;
  int i,j,k;
  int xSize, ySize, zSize;
  RealType cl, cp, cs;
  int cont=0;
  double def;
  double angle;

  pts = new vtkIdType[source->GetMaxCellSize()];
  trans = vtkTransform::New();
  matrix = vtkMatrix4x4::New();
//cout<<"3\n";
  // set up working matrices
  m[0] = m0; m[1] = m1; //m[2] = m2; 
  v[0] = v0; v[1] = v1; //v[2] = v2; 

//  vtkDebugMacro(<<"Generating tensor glyphs");

  pd = input->GetPointData();
  outPD = output->GetPointData();
//  inTensors = pd->GetTensors();
  inScalars = pd->GetScalars();
  numPts = input->GetNumberOfPoints();

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
  if ( newScalars )
    {
    int idx = outPD->AddArray(newScalars);
    outPD->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    newScalars->Delete();
    }

  if ( newNormals )
    {
    outPD->SetNormals(newNormals);
    newNormals->Delete();
    }


  //
  // Traverse all Input points, transforming glyph at Source points
  //
  trans->PreMultiply();

  for (inPtId=0; inPtId < numPts; inPtId++)
    {
    ptIncr = inPtId * numSourcePts;

    // Translation is postponed

//    inTensors->GetTuple(inPtId, tensor);
//cout<<"a\n";

    scalarValue = inScalars->GetTuple1(inPtId);

//ComputeEigenSystem(inPtId,v,w);

//if (tensor[0]!=0)
//cout<<v[0][0]<<" "<<v[1][0]<<"\t\t"<<v[0][1]<<" "<<v[1][1]<<"\t\t"<<w[0]<<" "<<w[1]<<" "<<"\n";


    //copy eigenvectors
//    xv[0] = v[0][0]; xv[1] = v[0][1];// xv[2] = v[2][0];
//    yv[0] = v[1][0]; yv[1] = v[1][1];// yv[2] = v[2][1];
//    zv[0] = v[0][2]; zv[1] = v[1][2]; zv[2] = v[2][2];

//if (tensor[0]!=0)
//cout<<xv[0]<<" "<<yv[0]<<"\t\t"<<xv[1]<<" "<<yv[1]<<"\t\t"<<w[0]<<" "<<w[1]<<" "<<"\n";

    // compute scale factors
//    w[0] *= this->ScaleFactor;
//    w[1] *= this->ScaleFactor;
//    w[2] *= this->ScaleFactor;


    if (scalarValue==0) continue;
/*
    if (fabs(w[0])>fabs(w[2]))
      factor = 1 / w[0];

    else factor = 1 / w[2];
*/
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


/*
cout<<pixel[0]<<" "<<pixel[1]<<" "<<pixel[2]<<" "<<pixel[3]<<" "<<pixel[4]<<" "<<pixel[5]<<"\n";
cout<<eigval[0]<<" "<<eigval[1]<<" "<<eigval[2]<<"\t";
cout<<eigvec[0][0]<<" "<<eigvec[0][1]<<" "<<eigvec[0][2]<<"\t";
cout<<eigvec[1][0]<<" "<<eigvec[1][1]<<" "<<eigvec[1][2]<<"\t";
cout<<eigvec[2][0]<<" "<<eigvec[2][1]<<" "<<eigvec[2][2]<<"\n";
cout<<w[0]<<" "<<w[1]<<" "<<w[2]<<"\t";
cout<<v[0][0]<<" "<<v[0][1]<<" "<<v[0][2]<<"\t";
cout<<v[1][0]<<" "<<v[1][1]<<" "<<v[1][2]<<"\t";
cout<<v[2][0]<<" "<<v[2][1]<<" "<<v[2][2]<<"\n\n";
*/

//    if ( (this->GlyphType==SUPERQUADRIC) || (this->ColorMode==COLOR_BY_CL) )
/*
    pixel.ComputeShapeCoefficients(cl,cp,cs);

    if (fabs(cs)>=csThreshold) continue;
*/
//if ( (w[0] < 0) || (w[1] < 0) || (w[2] < 0) )
//  cout<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<" "<<w[0]<<" "<<w[1]<<" "<<w[2]<<"\t"<<pixel[0]<<" "<<pixel[1]<<" "<<pixel[2]<<" "<<pixel[3]<<" "<<pixel[4]<<" "<<pixel[5]<<"\n";


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
    input->GetPoint(inPtId, x);
/*
    x[0] = origin[0] + spacing[0] * pixelIndex[0];
    x[1] = origin[1] + spacing[1] * pixelIndex[1];
    x[2] = origin[2] + spacing[2] * pixelIndex[2];
*/
    trans->Translate(x[0], x[1], x[2]);

    // normalized eigenvectors rotate object for eigen direction 0
 
//    double ang1 = atan2(yv[1],xv[1]) * 180 / 3.1415;

//cout<<"Angulo: "<<atan2(yv[0],xv[0])*180/3.1415<<" "<<atan2(yv[0],-xv[0])*180/3.1415<<"\n";
//cout<<"Angulo: "<<ang1<<"\n";

//    trans->RotateX(ang1);
//    trans->RotateY(ang1);

//    double signo = -1.0;
//    if (v[0][1]<0 && v[1][1]<0) signo = 1.0;
/*
    matrix->Element[0][0] = xv[1];
    matrix->Element[0][1] = 0;		//yv[0];
    matrix->Element[0][2] = yv[1];	//zv[0];
    matrix->Element[1][0] = 0;
    matrix->Element[1][1] = 1;
    matrix->Element[1][2] = 0;
    matrix->Element[2][0] = signo*xv[0];
    matrix->Element[2][1] = 0;		//yv[1];
    matrix->Element[2][2] = signo*yv[0];//zv[1];

    matrix->Element[0][0] = -1;
    matrix->Element[0][1] = 0;		//yv[0];
    matrix->Element[0][2] = -1;	//zv[0];
    matrix->Element[1][0] = 0;
    matrix->Element[1][1] = 1;
    matrix->Element[1][2] = 0;
    matrix->Element[2][0] = 1;
    matrix->Element[2][1] = 0;		//yv[1];
    matrix->Element[2][2] = -1;//zv[1];
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Ver página 64 del pdf
*/

    angle = atan2(-x[0],-x[2]);
    trans->RotateY(90 + angle*180/3.1415);


/* RECUPERAR ESTA PARTE SI SE SOLUCIONA EL CÁLCULO DE AUTOVECTORES
    matrix->Element[0][0] = eigvec[0][0];
    matrix->Element[0][1] = eigvec[1][0];
    matrix->Element[0][2] = eigvec[2][0];
    matrix->Element[1][0] = eigvec[0][1];
    matrix->Element[1][1] = eigvec[1][1];
    matrix->Element[1][2] = eigvec[2][1];
    matrix->Element[2][0] = eigvec[0][2];
    matrix->Element[2][1] = eigvec[1][2];
    matrix->Element[2][2] = eigvec[2][2];
*/
//    trans->Concatenate(matrix);

/* RECUPERAR ESTA PARTE SI SE SOLUCIONA EL CÁLCULO DE AUTOVECTORES
    trans->Scale(eigval[0], eigval[1], eigval[2]);
*/

//    trans->RotateY(90.0);

//    trans->Scale(factor*fabs(w[0]),factor*fabs(w[2]),.001);

//      trans->Scale(w[0],1,fabs(w[1]));
//if (fabs(w[0])>fabs(w[1])) factor = w[0];
//  else factor = w[1];
 
      trans->Scale(scalarValue,1,1);

    // multiply points (and normals if available) by resulting
    // matrix
    trans->TransformPoints(sourcePts,newPts);
    // Apply the transformation to a series of points, 
    // and append the results to outPts.
    if ( newNormals )
      {
      trans->TransformNormals(sourceNormals,newNormals);
      }

//if (fabs(w[0])>fabs(w[1])) scalarValue = w[0];
//  else scalarValue = w[1];


      // Copy point data from source
    for (i=0; i < numSourcePts; i++) 
      {
        newScalars->InsertNextTupleValue(&scalarValue);
      }
    cont++;
    }
//  vtkDebugMacro(<<"Generated " << numPts <<" tensor glyphs");

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

  trans->Delete();
  matrix->Delete();
//  source->Delete();

  return output;
}



vtkPolyData *vtkStrainTensorGlyph::GetSource()
{
  
  vtkPolyData *source;

  switch (this->GlyphType)
    {
    case ELLIPSOID: { vtkSphereSource *sphere = vtkSphereSource::New();
      sphere->SetPhiResolution(this->PhiResolution);
      sphere->SetThetaResolution(this->ThetaResolution);
      source = sphere->GetOutput();
      source->Update();
  cout<<"puntos: "<<source->GetNumberOfCells()<<"\n";
//      sphere->Delete();
      break;}

    case CUBOID: {vtkCubeSource *cube = vtkCubeSource::New();
      source = cube->GetOutput();
      source->Update();
//      cube->Delete();
      break;}

    case SUPERQUADRIC: {vtkSuperquadricSource *superquadric = vtkSuperquadricSource::New();
      superquadric->SetPhiResolution(this->PhiResolution);
      superquadric->SetThetaResolution(this->ThetaResolution);
      source = superquadric->GetOutput();
      source->Update();
//      superquadric->Delete();
      break;}
    }

  cout<<"puntos: "<<source->GetNumberOfCells()<<"\n";

  return source;

}

void vtkStrainTensorGlyph::ComputeEigenSystem (vtkIdType index, double *v[2], double w[2]) {

  double tensor[9];

  input->GetPointData()->GetTensors()->GetTuple(index, tensor);

  if(tensor[1]*tensor[2] <= 0.1e-20 ) {
    w[0] = tensor[0]; v[0][0] = 1; v[1][0] = 0;
    w[1] = tensor[3]; v[0][1] = 0; v[1][1] = 1;
    return;
  }

  // First of all, compute the eigenvalues, together with the rank of the filter:
  double tr = tensor[0] + tensor[3];
  double det = tensor[0] * tensor[3] - tensor[1] * tensor[2];
  double S = sqrt( (tr/2)*(tr/2) - det );
  w[0] = tr/2 + S;
  w[1] = tr/2 - S;
/*
  double SS = sqrt( std::max(((((tensor[0]-tensor[3])/2) + tensor[1] * tensor[2]) * ((tensor[0]-tensor[3])/2) + tensor[1] * tensor[2]), 0.0) );
  if( tensor[0] - tensor[3] < 0 ) {
    v[0][0] = tensor[1];
    v[1][0] = - (tensor[0]-tensor[3])/2 + SS;
    v[0][1] = + (tensor[0]-tensor[3])/2 - SS;
    v[1][1] = tensor[1];
  } else {
    v[0][1] = tensor[1];
    v[1][1] = - (tensor[0]-tensor[3])/2 - SS;
    v[0][0] = + (tensor[0]-tensor[3])/2 + SS;
    v[1][0] = tensor[1];
  }
*/

  if (tensor[1]!=0) {
    v[0][0]=w[0]-tensor[3];
    v[1][0]=tensor[1];
    v[0][1]=w[1]-tensor[3];
    v[1][1]=tensor[1];
  }

  double n1 = sqrt(v[0][0]*v[0][0]+v[1][0]*v[1][0]);
  v[0][0] /= n1; v[1][0] /= n1;

  double n2 = sqrt(v[0][1]*v[0][1]+v[1][1]*v[1][1]);
  v[0][1] /= n2; v[1][1] /= n2;

  return;

}









