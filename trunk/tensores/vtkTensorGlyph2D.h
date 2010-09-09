#ifndef __vtkTensorGlyph2D_h
#define __vtkTensorGlyph2D_h

#include "tensor/itkDTITensor.h"
#include <itkFixedArray.h>
#include <itkMatrix.h>
#include <itkNumericTraits.h>

#include "tensor/itkDTITensor.h"
#include "tensor/itkDWImages.h"

#include "strain/itkStrainTensor.h"

#include "vtkPolyData.h"

class vtkTensorGlyph2D
{

/*TODO
Parámetros de las fuentes (resolución phi y theta, ¿etc?)
Probar los set/get con macros

En .cxx
- Definir y crear la fuente
- escalares - hay q calcularlos (inScalars¿¿??)
- Situar los glifos
*/

public:

  enum { Dimension =  4 };

  typedef itk::NumericTraits<float>::RealType 			RealType;
  typedef itk::StrainTensor<float>				StrainPixelType;
  typedef itk::Image<StrainPixelType, Dimension>		StrainImageType;
//  typedef typename itk::NumericTraits<TComponent>::RealType RealType;
  typedef itk::FixedArray<RealType,2>                       EigenValuesArrayType;
  typedef itk::Matrix<RealType,2,2>                         EigenVectorsMatrixType;
  typedef itk::FixedArray<RealType,2>                       EigenVectorType;

  typedef enum ColorModes{INV,EIG0,EIG1,ST0,ST1,ST2} ColorModes;

  static vtkTensorGlyph2D *New();
  vtkTensorGlyph2D();
  ~vtkTensorGlyph2D();


  void SetInput(StrainImageType::Pointer entrada)
    {this->input = entrada;};

  vtkPolyData *GetOutput();
  vtkPolyData *GetDeformOutput();

  void SetGlyphType(int type)
    {this->GlyphType = type;};
  int GetGlyphType()
    {return this->GlyphType;};

  void SetSource(vtkPolyData *src) 
    {this->source = src;};
  vtkPolyData *GetSource();

  void SetScaling(int value)
    {this->Scaling = value;};
  int GetScaling()
    {return this->Scaling;};

  void SetScaleFactor(double value)
    {this->ScaleFactor = value;};
  double GetScaleFactor()
    {return this->ScaleFactor;};

  void SetColorMode(int value)
    {this->ColorMode = value;};
  int GetColorMode()
    {return this->ColorMode;};

  void SetPhiResolution(int value)
    {this->PhiResolution = value;};
  int GetPhiResolution()
    {return this->PhiResolution;};

  void SetThetaResolution(int value)
    {this->ThetaResolution = value;};
  int GetThetaResolution()
    {return this->ThetaResolution;};

  void SetcsThreshold(double value)
    {this->csThreshold = value;};
  double GetcsThreshold()
    {return this->csThreshold;};

  void SetPlanoZ(int value)
    {this->PlanoZ = value;};
  int GetPlanoZ()
    {return this->PlanoZ;};

  void SetTiempo(int value)
    {this->Tiempo = value;};
  int GetTiempo()
    {return this->Tiempo;};


  void ComputeEigenSystem(vtkIdType,double*[2],double[2]);

//  void SetBounds(int *values)
//    {this->Bounds = values;};
  void SetBounds(int v1,int v2,int v3,int v4,int v5,int v6)
  {
    Bounds[0]=v1;
    Bounds[1]=v2;
    Bounds[2]=v3;
    Bounds[3]=v4;
    Bounds[4]=v5;
    Bounds[5]=v6;
  };
  int *GetBounds()
    {return this->Bounds;};

protected:

  StrainImageType::Pointer input;
  vtkPolyData *source;
  vtkPoints *inputPoints;
  int GlyphType;

  int Scaling; // Determine whether scaling of geometry is performed
  double ScaleFactor; // Scale factor to use to scale geometry
  int ColorMode; // The coloring mode to use for the glyphs.

  int PlanoZ, Tiempo;

  int Bounds[6];

  int PhiResolution;
  int ThetaResolution;

  double csThreshold;

};

#endif
