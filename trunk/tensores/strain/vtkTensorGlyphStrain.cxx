/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkTensorGlyphStrain.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTensorGlyphStrain.h"

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

//vtkCxxRevisionMacro(vtkTensorGlyphStrain, "$Revision: 1.57.12.1 $");
//vtkStandardNewMacro(vtkTensorGlyphStrain);

// Construct object with scaling on and scale factor 1.0. Eigenvalues are 
// extracted, glyphs are colored with input scalar data, and logarithmic
// scaling is turned off.
vtkTensorGlyphStrain::vtkTensorGlyphStrain()
{
  this->inputPoints = NULL;

  this->Scaling = 1;
  this->ScaleFactor = 1.0;
  this->ColorMode = DEF;

}

vtkTensorGlyphStrain::~vtkTensorGlyphStrain()
{
}

vtkPolyData *vtkTensorGlyphStrain::GetOutput()
{

  vtkPolyData *source;
  vtkPolyData *output = vtkPolyData::New();

  vtkIdType numPts, numSourcePts, numSourceCells, inPtId;
  vtkPoints *sourcePts;
  vtkDataArray *sourceNormals;
  vtkCellArray *sourceCells, *cells;  
  vtkPoints *newPts;
  vtkFloatArray *newScalars=NULL;
  vtkFloatArray *newNormals=NULL;
  double x[3];
  vtkTransform *trans;
  vtkCell *cell;
  vtkIdList *cellPts;
  int npts;
  vtkIdType *pts;
  vtkIdType ptIncr, cellId;
  vtkMatrix4x4 *matrix;
  double maxScale;
  vtkPointData *pd, *outPD;

  STImageType::IndexType pixelIndex;
  STPixelType pixel;
  EigenValuesArrayType eigval;
  EigenVectorsMatrixType eigvec;
  double factor;
  float scalarValue;
  int i,j,k;
  int xSize, ySize;
  int cont=0;
  double angulo;

  double tensor1[3],tensor2[3],tensor3[3],tensor4[3];

  RealType t1,t2;

  STImageType::PointType origin = this->input->GetOrigin();
  STImageType::SpacingType spacing = this->input->GetSpacing();

  vtkGlyphSource2D *glyphsource = vtkGlyphSource2D::New();
  glyphsource->SetGlyphTypeToDiamond();
  glyphsource->FilledOff();
  glyphsource->Update();
  source = glyphsource->GetOutput();
  source->Update();

  pts = new vtkIdType[source->GetMaxCellSize()];
  trans = vtkTransform::New();
  matrix = vtkMatrix4x4::New();

  outPD = output->GetPointData();

  xSize = this->input->GetRequestedRegion().GetSize()[0];
  ySize = this->input->GetRequestedRegion().GetSize()[1];
  numPts = xSize * ySize;

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

//double factor = 0;
factor = 0;
double deform_max = 0;
  for (inPtId=0; inPtId<numPts; inPtId++)
    {
    ptIncr = inPtId * numSourcePts;

    // Translation is postponed
    pixelIndex[0] = inPtId % xSize;
    pixelIndex[1] = inPtId / xSize;
    pixelIndex[2] = this->PlanoZ;
    pixelIndex[3] = this->Tiempo;
    pixel = input->GetPixel(pixelIndex);

    pixel.ComputeEigenSystem(eigval,eigvec);

    if (eigval[0]==0) continue;
    if (fabs(eigval[0])>factor) factor = fabs(eigval[0]);
    if (fabs(eigval[1])>factor) factor = fabs(eigval[1]);

    int indice = inPtId + PlanoZ*56*56*16 + Tiempo*56*56;
//    double deform = sqrt (deformArray->GetTuple(indice)[0]*deformArray->GetTuple(indice)[0] + deformArray->GetTuple(indice)[1]*deformArray->GetTuple(indice)[1] );
    double deform = sqrt (deformImage->GetPixel(pixelIndex)[0]*deformImage->GetPixel(pixelIndex)[0] + deformImage->GetPixel(pixelIndex)[1]*deformImage->GetPixel(pixelIndex)[1] );
    if (deform>deform_max) deform_max = deform;
  }

factor = 1/factor;

  //
  // Traverse all Input points, transforming glyph at Source points
  //
  trans->PreMultiply();

  for (inPtId=0; inPtId<numPts; inPtId++)
    {
    ptIncr = inPtId * numSourcePts;

    // Translation is postponed
    pixelIndex[0] = inPtId % xSize;
    pixelIndex[1] = inPtId / xSize;
    pixelIndex[2] = this->PlanoZ;
    pixelIndex[3] = this->Tiempo;
    pixel = input->GetPixel(pixelIndex);

    pixel.ComputeEigenSystem(eigval,eigvec);

    if (eigval[0]==0) continue;

    // Now do the real work for each "direction"

    // Remove previous scales ...
    trans->Identity();

    // translate Source to Input point
    x[0] = origin[0] + spacing[0] * pixelIndex[0];
    x[1] = origin[1] + spacing[1] * pixelIndex[1];
    x[2] = origin[2] + spacing[2] * pixelIndex[2];
//    trans->Translate(x[0], x[1], x[2]);
    trans->Translate(x[0]+spacing[0]/2, x[1]+spacing[1]/2, x[2]);

    angulo = atan(eigvec[1][0]/eigvec[0][0]) * 180 / 3.1415;
    trans->RotateZ(angulo);

//    if (fabs(eigval[0])>fabs(eigval[1])) factor = fabs(1 / eigval[0]);
//      else factor = fabs(1 / eigval[1]);

//    double factor1 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[0]*factor+1) + 0.5*sqrt(eigval[0]*eigval[0]*factor*factor+eigval[1]*eigval[1]*factor*factor));
//    double factor2 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[1]*factor+1) + 0.5*sqrt(eigval[0]*eigval[0]*factor*factor+eigval[1]*eigval[1]*factor*factor));

//    int indice = inPtId + PlanoZ*56*56*16 + Tiempo*56*56;
//    double deform = sqrt (deformArray->GetTuple(indice)[0]*deformArray->GetTuple(indice)[0] + deformArray->GetTuple(indice)[1]*deformArray->GetTuple(indice)[1] );

//if (inPtId == 11954) {
//cout<<pixelIndex<<" deform1 "<<deform<<" deform_max "<<deform_max<<"\n";
//cout<<"deform viejo "<<deformArray->GetTuple(indice)[0]<<" deform_max "<<deformArray->GetTuple(indice)[1]<<"\n";
//cout<<"deform nuevo "<<deformImage->GetPixel(pixelIndex)[0]<<" deform_max "<<deformImage->GetPixel(pixelIndex)[1]<<"\n";
//}
//cout<<"deform final "<<deform<<"\n";

    double deform = sqrt (deformImage->GetPixel(pixelIndex)[0]*deformImage->GetPixel(pixelIndex)[0] + deformImage->GetPixel(pixelIndex)[1]*deformImage->GetPixel(pixelIndex)[1] );
    deform = deform / deform_max;

    double factor1 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[0]*factor+1) + 0.5*deform);
    double factor2 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[1]*factor+1) + 0.5*deform);

    trans->Scale(this->ScaleFactor*factor1,this->ScaleFactor*factor2,1);

    // multiply points (and normals if available) by resulting matrix
    trans->TransformPoints(sourcePts,newPts);

    // Apply the transformation to a series of points, 
    // and append the results to outPts.
    if ( newNormals )
      {
      trans->TransformNormals(sourceNormals,newNormals);
      }

//    if (pixelIndex[0]==28) cout<<pixelIndex<<" pixel "<<pixel<<" eigval "<<eigval<<" invariante "<<pixel.GetInvariant()<<" deform "<<deformImage->GetPixel(pixelIndex)[0]<<" "<<deformImage->GetPixel(pixelIndex)[1]<<" "<<deform<<" max_def "<<deform_max<<"\n";

//    if (fabs(eigval[0])>fabs(eigval[1])) scalarValue = eigval[0];
//      else scalarValue = eigval[1];

    switch (this->ColorMode)
    {
        case DEF: 
          scalarValue = deform * deform_max;
//          scalarValue = pixel.GetInvariant();
          break;

        case EIG0: 
          scalarValue = eigval[0];
          break;

        case EIG1: 
          scalarValue = eigval[1];
          break;

        case ST0: 
          scalarValue = pixel[0];
          break;

        case ST1: 
          scalarValue = pixel[1];
          break;

        case ST2: 
          scalarValue = pixel[2];
          break;

    }

      // Copy point data from source
    for (i=0; i < numSourcePts; i++) 
      {
        newScalars->InsertNextTupleValue(&scalarValue);
      }

    cont++;

    }







  if (this->inputPoints) {
  
  numPts = this->inputPoints->GetNumberOfPoints();

  for (inPtId=0; inPtId<numPts; inPtId++)
  {
    ptIncr = inPtId * numSourcePts;

    this->inputPoints->GetPoint(inPtId,x);

    pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
    pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);
    pixel = input->GetPixel(pixelIndex);

    double deform_values[2];
    double eigvalues[2];
    interpolarTensor (x, eigvalues, &angulo, deform_values);

    if (eigvalues[0] == 0) continue;
     // Now do the real work for each "direction"

    // Remove previous scales ...
    trans->Identity();

//    trans->Translate(x[0], x[1], x[2]);
    trans->Translate(x[0]+spacing[0]/2, x[1]+spacing[1]/2, x[2]);

    trans->RotateZ(angulo);

//    if (fabs(eigval[0])>fabs(eigval[1])) factor = fabs(1 / eigval[0]);
//      else factor = fabs(1 / eigval[1]);

//    eigval[0] *= factor;
//    eigval[1] *= factor;

//    double factor1 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[0]*factor+1) + 0.5*sqrt(eigval[0]*eigval[0]*factor*factor+eigval[1]*eigval[1]*factor*factor));
//    double factor2 = sqrt(2)*(0.25 + 0.25*0.5*(eigval[1]*factor+1) + 0.5*sqrt(eigval[0]*eigval[0]*factor*factor+eigval[1]*eigval[1]*factor*factor));

    double deform = sqrt (deform_values[0]*deform_values[0] + deform_values[1]*deform_values[1] );
    deform = deform / deform_max;

    double factor1 = sqrt(2)*(0.25 + 0.25*0.5*(eigvalues[0]*factor+1) + 0.5*deform);
    double factor2 = sqrt(2)*(0.25 + 0.25*0.5*(eigvalues[1]*factor+1) + 0.5*deform);

    trans->Scale(this->ScaleFactor*factor1,this->ScaleFactor*factor2,1);
//  trans->Scale(this->ScaleFactor*factor*eigval[0],this->ScaleFactor*factor*eigval[1],1);

    // multiply points (and normals if available) by resulting matrix
    trans->TransformPoints(sourcePts,newPts);

    // Apply the transformation to a series of points, 
    // and append the results to outPts.
    if ( newNormals )
      {
      trans->TransformNormals(sourceNormals,newNormals);
      }

//    if (fabs(eigval[0])>fabs(eigval[1])) scalarValue = eigval[0];
//      else scalarValue = eigval[1];

    switch (this->ColorMode)
    {
        case DEF: 
          scalarValue = deform * deform_max;
//          scalarValue = eigvalues[0]*eigvalues[0] + eigvalues[1]*eigvalues[1]; //pixel.GetInvariant();
          break;

        case EIG0: 
          scalarValue = eigvalues[0];
          break;

        case EIG1: 
          scalarValue = eigvalues[1];
          break;

        case ST0: 
          scalarValue = pixel[0];
          break;

        case ST1: 
          scalarValue = pixel[1];
          break;

        case ST2: 
          scalarValue = pixel[2];
          break;

    }

      // Copy point data from source
    for (i=0; i < numSourcePts; i++) 
      {
        newScalars->InsertNextTupleValue(&scalarValue);
      }

    cont++;
    }
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
      for (i=0; i < npts; i++)
        {
        pts[i] = cellPts->GetId(i) + ptIncr;
        }
      output->InsertNextCell(cell->GetCellType(),npts,pts);
      }
    }

  //
  // Update output and release memory
  //
  delete [] pts;

  output->SetPoints(newPts);

  glyphsource->Delete();
  newPts->Delete();
  trans->Delete();
  matrix->Delete();

  return output;
}

void vtkTensorGlyphStrain::interpolarTensor (double x[3], double eigval_out[2], double *angulo, double deform_values[2]) {

    double eigval0_1, eigval0_2, eigval0_3, eigval0_4;
    double eigval1_1, eigval1_2, eigval1_3, eigval1_4;
    double angulo1, angulo2, angulo3, angulo4;
    double deform0_1, deform0_2, deform0_3, deform0_4;
    double deform1_1, deform1_2, deform1_3, deform1_4;

    double t1, t2;

    STPixelType pixel;
    DeformPixelType deformPixel;

    STImageType::PointType origin = this->input->GetOrigin();
    STImageType::SpacingType spacing = this->input->GetSpacing();
    STImageType::IndexType pixelIndex;

    EigenValuesArrayType eigval;
    EigenVectorsMatrixType eigvec;

    pixelIndex[0] = (int) ((x[0]-origin[0]) / spacing[0]);
    pixelIndex[1] = (int) ((x[1]-origin[1]) / spacing[1]);

    t1 = ( (x[0]-origin[0]) / spacing[0] - pixelIndex[0] );// / spacing[0];
    t2 = ( (x[1]-origin[1]) / spacing[1] - pixelIndex[1] );// / spacing[1];

//    pixelIndex[0] = round( (x[0]-origin[0]) / spacing[0] ); //inPtId % xSize;
//    pixelIndex[1] = round( (x[1]-origin[1]) / spacing[1] ); //inPtId / xSize;
    pixelIndex[2] = this->PlanoZ;
    pixelIndex[3] = this->Tiempo;
    pixel = input->GetPixel(pixelIndex);
    deformPixel = deformImage->GetPixel(pixelIndex);
    deform0_1 = deformPixel[0]; deform1_1 = deformPixel[1];
    pixel.ComputeEigenSystem(eigval,eigvec);
    eigval0_1=eigval[0]; eigval1_1=eigval[1]; angulo1=atan(eigvec[1][0]/eigvec[0][0]);
    double angulo0 = angulo1 - (3.1415/2);
    angulo1 -= angulo0;
    deformPixel = deformImage->GetPixel(pixelIndex);
    deform0_1 = deformPixel[0]; deform1_1 = deformPixel[1];

    pixelIndex[0] += 1;
    pixel = input->GetPixel(pixelIndex);
    pixel.ComputeEigenSystem(eigval,eigvec);
    eigval0_2=eigval[0]; eigval1_2=eigval[1]; angulo2=atan(eigvec[1][0]/eigvec[0][0]);
    angulo2 -= angulo0;
    if (angulo2 > 3.1415) angulo2 -= 3.1415;
    if (angulo2 < 0) angulo2 += 3.1415;
    deformPixel = deformImage->GetPixel(pixelIndex);
    deform0_2 = deformPixel[0]; deform1_2 = deformPixel[1];
	
    pixelIndex[0] -= 1;
    pixelIndex[1] += 1;
    pixel = input->GetPixel(pixelIndex);
    pixel.ComputeEigenSystem(eigval,eigvec);
    eigval0_3=eigval[0]; eigval1_3=eigval[1]; angulo3=atan(eigvec[1][0]/eigvec[0][0]);
    angulo3 -= angulo0;
    if (angulo3 > 3.1415) angulo3 -= 3.1415;
    if (angulo3 < 0) angulo3 += 3.1415;
    deformPixel = deformImage->GetPixel(pixelIndex);
    deform0_3 = deformPixel[0]; deform1_3 = deformPixel[1];

    pixelIndex[0] += 1;
    pixel = input->GetPixel(pixelIndex);
    pixel.ComputeEigenSystem(eigval,eigvec);
    eigval0_4=eigval[0]; eigval1_4=eigval[1]; angulo4=atan(eigvec[1][0]/eigvec[0][0]);
    angulo4 -= angulo0;
    if (angulo4 > 3.1415) angulo4 -= 3.1415;
    if (angulo4 < 0) angulo4 += 3.1415;
    deformPixel = deformImage->GetPixel(pixelIndex);
    deform0_4 = deformPixel[0]; deform1_4 = deformPixel[1];

// VER QUÉ PASA SI ALGUNO ES CERO
//    if ( (tensor1[0]==0) || (tensor2[0]==0) || (tensor3[0]==0) || (tensor4[0]==0) ) continue;
//    if ( (tensor1[1]==0) || (tensor2[1]==0) || (tensor3[1]==0) || (tensor4[1]==0) ) continue;

    eigval_out[0] = 0;
    eigval_out[1] = 0;
    *angulo = 0;
    deform_values[0] = 0;
    deform_values[1] = 0;

    if (t1 == 0 && t2 == 0) {
       eigval_out[0] = eigval0_1;
       eigval_out[1] = eigval1_1;

       deform_values[0] = deform0_1;
       deform_values[1] = deform1_1;

       *angulo = angulo1;
    }

    else if (t1 == 0) {
      if (eigval0_1 == 0 || eigval0_3 == 0) return;

      eigval_out[0] = (1-t2) * eigval0_1 + t2 * eigval0_3;
      eigval_out[1] = (1-t2) * eigval1_1 + t2 * eigval1_3;

      deform_values[0] = (1-t2) * deform0_1 + t2 * deform0_3;
      deform_values[1] = (1-t2) * deform1_1 + t2 * deform1_3;

      *angulo = (1-t2) * angulo1 + t2 * angulo3;
    } 

    else if (t2 ==0) {
      if (eigval0_1 == 0 || eigval0_2 == 0) return;

      eigval_out[0] = (1-t1) * eigval0_1 + t1 * eigval0_2;
      eigval_out[1] = (1-t1) * eigval1_1 + t1 * eigval1_2;

      deform_values[0] = (1-t1) * deform0_1 + t1 * deform0_2;
      deform_values[1] = (1-t1) * deform1_1 + t1 * deform1_2;

      *angulo = (1-t1) * angulo1 + t1 * angulo2;
    } 

    else if (t1 != 0 && t2 != 0) {
      if (eigval0_1 == 0 || eigval0_2 == 0 || eigval0_3 == 0 || eigval0_4 == 0) return;

      eigval_out[0] = (1-t2) * ((1-t1)*eigval0_1 + t1*eigval0_2) + t2 * ((1-t1)*eigval0_3 + t1*eigval0_4);
      eigval_out[1] = (1-t2) * ((1-t1)*eigval1_1 + t1*eigval1_2) + t2 * ((1-t1)*eigval1_3 + t1*eigval1_4);

      deform_values[0] = (1-t2) * ((1-t1)*deform0_1 + t1*deform0_2) + t2 * ((1-t1)*deform0_3 + t1*deform0_4);
      deform_values[1] = (1-t2) * ((1-t1)*deform1_1 + t1*deform1_2) + t2 * ((1-t1)*deform1_3 + t1*deform1_4);

      *angulo = (1-t2) * ((1-t1)*angulo1 + t1*angulo2) + t2 * ((1-t1)*angulo3 + t1*angulo4);
    }

    *angulo = (*angulo + angulo0) * 180 / 3.1415;


}


