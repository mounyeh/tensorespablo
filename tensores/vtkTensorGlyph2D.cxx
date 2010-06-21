/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkTensorGlyph2D.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTensorGlyph2D.h"

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
#include "vtkGlyphSource2D.h"

#include <algorithm>

//vtkCxxRevisionMacro(vtkTensorGlyph2D, "$Revision: 1.57.12.1 $");
//vtkStandardNewMacro(vtkTensorGlyph2D);

// Construct object with scaling on and scale factor 1.0. Eigenvalues are 
// extracted, glyphs are colored with input scalar data, and logarithmic
// scaling is turned off.
vtkTensorGlyph2D::vtkTensorGlyph2D()
{
  this->inputPoints = NULL;

  this->Scaling = 1;
  this->ScaleFactor = 1.0;
  this->ColorMode = COLOR_BY_SCALARS;
  this->ClampScaling = 0;
  this->MaxScaleFactor = 100;

  this->PhiResolution = 8;
  this->ThetaResolution = 8;

  this->Bounds[0] = 0;
  this->Bounds[1] = 0;
  this->Bounds[2] = 0;
  this->Bounds[3] = 0;
  this->Bounds[4] = 0;
  this->Bounds[5] = 0;

  this->csThreshold = 1;

  vtkGlyphSource2D *glyphsource = vtkGlyphSource2D::New();
  glyphsource->SetGlyphTypeToDiamond();
  glyphsource->SetScale(100);
  glyphsource->FilledOff();
  glyphsource->Update();
  this->source = glyphsource->GetOutput();
//  cout<<this->source->GetNumberOfCells()<<"\n";
//  glyphsource->Delete();
//  cout<<this->source->GetNumberOfCells()<<"\n";
}

vtkTensorGlyph2D::~vtkTensorGlyph2D()
{
}

vtkPolyData *vtkTensorGlyph2D::GetOutput()
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
//  input->Update();
//cout<<"2\n";
  vtkDataArray *inTensors;
  double tensor[4];
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
  double xv[2], yv[2];
  double maxScale;
  vtkPointData *pd, *outPD;

  StrainImageType::IndexType pixelIndex;
  StrainPixelType pixel;
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
  double angulo;

  double tensor1[4],tensor2[4],tensor3[4],tensor4[4];
  double info1[3],info2[3],info3[3],info4[3];

  StrainImageType::PointType origin = this->input->GetOrigin();
  StrainImageType::SpacingType spacing = this->input->GetSpacing();

cout<<"Origen: "<<origin<<"\tEspaciado: "<<spacing<<"\n";

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

//  pd = input->GetPointData();
  outPD = output->GetPointData();
//  inTensors = pd->GetTensors();
//  inScalars = pd->GetScalars();
//  numPts = input->GetNumberOfPoints();

  xSize = this->input->GetRequestedRegion().GetSize()[0];
  ySize = this->input->GetRequestedRegion().GetSize()[1];
  numPts = xSize * ySize;

  source->Update();

  //
  // Allocate storage for output PolyData
  //

  sourcePts = source->GetPoints();
  numSourcePts = sourcePts->GetNumberOfPoints();
  numSourceCells = source->GetNumberOfCells();
//cout<<"5\n";
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
//int indiceInicial = 100 * this->PlanoZ + 10 * this->Tiempo;
//cout<<"8\n";
//  for (inPtId=0; inPtId < numPts; inPtId++)
  for (inPtId=0; inPtId<numPts; inPtId++)
    {
    ptIncr = inPtId * numSourcePts;

    // Translation is postponed

    pixelIndex[0] = inPtId % xSize;
    pixelIndex[1] = inPtId / xSize;
    pixelIndex[2] = this->PlanoZ;
    pixelIndex[3] = this->Tiempo;
    pixel = input->GetPixel(pixelIndex);

//    inTensors->GetTuple(inPtId, tensor);
//cout<<"a\n";
//cout<<v[1][0]<<"aa\n";
//    scalarValue = inScalars->GetTuple1(inPtId);

//cout<<v[1][0]<<"cc\n";
    pixel.ComputeEigenSystem(eigval,eigvec);

//cout<<"c\n";
//cout<<tensor[0]<<" "<<tensor[1]<<" "<<tensor[2]<<" "<<tensor[3]<<"\n";
//cout<<v[0][0]<<" "<<v[1][0]<<"\t\t"<<v[0][1]<<" "<<v[1][1]<<"\t\t"<<w[0]<<" "<<w[1]<<" "<<"\n";


    //copy eigenvectors
//    xv[0] = v[0][0]; xv[1] = v[0][1];// xv[2] = v[2][0];
//    yv[0] = v[1][0]; yv[1] = v[1][1];// yv[2] = v[2][1];
//    zv[0] = v[0][2]; zv[1] = v[1][2]; zv[2] = v[2][2];

    // compute scale factors
//    eigval[0] *= this->ScaleFactor;
//    eigval[1] *= this->ScaleFactor;

//    if (scalarValue==0) continue;

//cout<<"e\n";

     // Now do the real work for each "direction"

    // Remove previous scales ...
    trans->Identity();

    // translate Source to Input point
//    input->GetPoint(inPtId, x);
//    v[0]=v0; v[1]=v1;

    x[0] = origin[0] + spacing[0] * pixelIndex[0];
    x[1] = origin[1] + spacing[1] * pixelIndex[1];
    x[2] = origin[2] + spacing[2] * pixelIndex[2];

    trans->Translate(x[0], x[1], x[2]);
//cout<<"Puntos: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<"\n";

    // normalized eigenvectors rotate object for eigen direction 0
 
    angulo = atan(eigvec[1][0]/eigvec[0][0]) * 180 / 3.1415;
    trans->RotateZ(angulo);

//cout<<"g\n";

    if (fabs(eigval[0])>fabs(eigval[1])) factor = fabs(1 / eigval[0]);
      else factor = fabs(1 / eigval[1]);

    eigval[0] *= factor;
    eigval[1] *= factor;

    double factor1 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[0]+1) + 0.5*sqrt(eigval[0]*eigval[0]+eigval[1]*eigval[1]));
    double factor2 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[1]+1) + 0.5*sqrt(eigval[0]*eigval[0]+eigval[1]*eigval[1]));
//    double factor1 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[0]+1));
//    double factor2 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[1]+1));

if (inPtId==0) cout<<"Factores: "<<factor1<<" "<<factor2<<"\n";

    trans->Scale(this->ScaleFactor*factor1,this->ScaleFactor*factor2,1);
//    trans->Scale(0.7,0.7,1);

    // multiply points (and normals if available) by resulting
    // matrix
    trans->TransformPoints(sourcePts,newPts);

    // Apply the transformation to a series of points, 
    // and append the results to outPts.
    if ( newNormals )
      {
      trans->TransformNormals(sourceNormals,newNormals);
      }

    if (fabs(eigval[0])>fabs(eigval[1])) scalarValue = eigval[0];
      else scalarValue = eigval[1];

    scalarValue = 1;

      // Copy point data from source
    for (i=0; i < numSourcePts; i++) 
      {
        newScalars->InsertNextTupleValue(&scalarValue);
      }
    cont++;
    }
//  vtkDebugMacro(<<"Generated " << numPts <<" tensor glyphs");








// ** Descomentar estas líneas para permitir puntos de entrada **
//  if (this->inputPoints) numPts = this->inputPoints->GetNumberOfPoints();
//  if (this->inputPoints)


  for (inPtId=0; inPtId<numPts; inPtId++)
    {
    ptIncr = inPtId * numSourcePts;

// ** Descomentar estas líneas para permitir puntos de entrada **

//    this->inputPoints->GetPoint(inPtId,x);

//    pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
//    pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);

//    t1 = ( (x[0]-origin[0]) / spacing[0] - pixelIndex[0] ) / spacing[0];
//    t2 = ( (x[1]-origin[1]) / spacing[1] - pixelIndex[1] ) / spacing[1];



    pixelIndex[0] = inPtId % xSize;
    pixelIndex[1] = inPtId / xSize;
    pixelIndex[2] = this->PlanoZ;
    pixelIndex[3] = this->Tiempo;
    pixel = input->GetPixel(pixelIndex);

    pixel.ComputeEigenSystem(eigval,eigvec);
    info1[0]=eigval[0]; info1[1]=eigval[1]; info1[2]=atan(eigvec[1][0]/eigvec[0][0]);
// OJO!! A VER SI CALCULO BIEN EL ANGULO
//info1[2]=atan2(yv[1],xv[1]);

//cout<<"a "<<eigvec[0][0]<<" "<<eigvec[1][0]<<"\t\t"<<eigvec[0][1]<<" "<<eigvec[1][1]<<"\t\t"<<eigval[0]<<" "<<eigval[1]<<" "<<"\n";

    pixelIndex[0] += 1;
    pixel = input->GetPixel(pixelIndex);
    pixel.ComputeEigenSystem(eigval,eigvec);
    info2[0]=eigval[0]; info2[1]=eigval[1]; info2[2]=atan(eigvec[1][0]/eigvec[0][0]);

    pixelIndex[0] -= 1;
    pixelIndex[1] += 1;
    pixel = input->GetPixel(pixelIndex);
    pixel.ComputeEigenSystem(eigval,eigvec);
    info3[0]=eigval[0]; info3[1]=eigval[1]; info3[2]=atan(eigvec[1][0]/eigvec[0][0]);

    pixelIndex[0] += 1;
    pixel = input->GetPixel(pixelIndex);
    pixel.ComputeEigenSystem(eigval,eigvec);
    info4[0]=eigval[0]; info4[1]=eigval[1]; info4[2]=atan(eigvec[1][0]/eigvec[0][0]);

// ** Descomentar estas líneas para permitir puntos de entrada **
//    eigval[0] = t2 * (t1*info1[0] + (1-t1)*info2[0]) + (1-t2) * (t1*info3[0] + (1-t1)*info4[0]);
//    eigval[1] = t2 * (t1*info1[1] + (1-t1)*info2[1]) + (1-t2) * (t1*info3[1] + (1-t1)*info4[1]);
//    angulo = ( t2 * (t1*info1[0] + (1-t1)*info2[0]) + (1-t2) * (t1*info3[0] + (1-t1)*info4[0]) ) * 180 / 3.1415;

    eigval[0] = 0.25 * (info1[0]+info2[0]+info3[0]+info4[0]);
    eigval[1] = 0.25 * (info1[1]+info2[1]+info3[1]+info4[1]);
    angulo = 0.25 * (info1[2]+info2[2]+info3[2]+info4[2]) * 180 / 3.1415;


//    inTensors->GetTuple(inPtId, tensor);
//cout<<"a\n";
//cout<<v[1][0]<<"aa\n";
//    scalarValue = inScalars->GetTuple1(inPtId);

//cout<<v[1][0]<<"cc\n";
//    pixel.ComputeEigenSystem(eigval,eigvec);

//cout<<"c\n";
//cout<<tensor[0]<<" "<<tensor[1]<<" "<<tensor[2]<<" "<<tensor[3]<<"\n";
//cout<<"b "<<eigvec[0][0]<<" "<<eigvec[1][0]<<"\t\t"<<eigvec[0][1]<<" "<<eigvec[1][1]<<"\t\t"<<eigval[0]<<" "<<eigval[1]<<" "<<"\n";


    //copy eigenvectors
//    xv[0] = v[0][0]; xv[1] = v[0][1];// xv[2] = v[2][0];
//    yv[0] = v[1][0]; yv[1] = v[1][1];// yv[2] = v[2][1];
//    zv[0] = v[0][2]; zv[1] = v[1][2]; zv[2] = v[2][2];

    // compute scale factors
//    eigval[0] *= this->ScaleFactor;
//    eigval[1] *= this->ScaleFactor;

//    if (scalarValue==0) continue;

    if (fabs(eigval[0])>fabs(eigval[2]))
      factor = 1 / eigval[0];

    else factor = 1 / eigval[2];
//cout<<"e\n";

     // Now do the real work for each "direction"

    // Remove previous scales ...
    trans->Identity();

    // translate Source to Input point
//    input->GetPoint(inPtId, x);
//    v[0]=v0; v[1]=v1;

    x[0] = origin[0] + spacing[0] * pixelIndex[0] - spacing[0]/2;
    x[1] = origin[1] + spacing[1] * pixelIndex[1] - spacing[1]/2;
    x[2] = origin[2] + spacing[2] * pixelIndex[2];

    trans->Translate(x[0], x[1], x[2]);
//cout<<"Puntos: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<"\n";

    // normalized eigenvectors rotate object for eigen direction 0
 
//    angulo = atan2(eigvec[1][1],eigvec[0][1]) * 180 / 3.1415;
    trans->RotateZ(angulo);

//cout<<"g\n";

    if (fabs(eigval[0])>fabs(eigval[1])) factor = fabs(1 / eigval[0]);
      else factor = fabs(1 / eigval[1]);

    eigval[0] *= factor;
    eigval[1] *= factor;

    double factor1 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[0]+1) + 0.5*sqrt(eigval[0]*eigval[0]+eigval[1]*eigval[1]));
    double factor2 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[1]+1) + 0.5*sqrt(eigval[0]*eigval[0]+eigval[1]*eigval[1]));
//    double factor1 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[0]+1));
//    double factor2 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[1]+1));

    trans->Scale(this->ScaleFactor*factor1,this->ScaleFactor*factor2,1);
//    trans->Scale(0.7,0.7,1);

    // multiply points (and normals if available) by resulting
    // matrix
    trans->TransformPoints(sourcePts,newPts);

    // Apply the transformation to a series of points, 
    // and append the results to outPts.
    if ( newNormals )
      {
      trans->TransformNormals(sourceNormals,newNormals);
      }

    if (fabs(eigval[0])>fabs(eigval[1])) scalarValue = eigval[0];
      else scalarValue = eigval[1];

    scalarValue = .2;

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

void vtkTensorGlyph2D::ComputeEigenSystem (vtkIdType index, double *v[2], double w[2]) {

  double tensor[9];

//  input->GetPointData()->GetTensors()->GetTuple(index, tensor);

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


