#ifndef __vtkTensorGlyphStrain_h
#define __vtkTensorGlyphStrain_h

#include "tensor/itkDTITensor.h"
#include <itkFixedArray.h>
#include <itkMatrix.h>
#include <itkNumericTraits.h>

#include "tensor/itkDTITensor.h"
#include "tensor/itkDWImages.h"

#include "strain/itkStrainTensor.h"

#include "vtkPolyData.h"
#include "vtkFloatArray.h"

class vtkTensorGlyphStrain
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

  typedef itk::NumericTraits<float>::RealType 		RealType;
  typedef itk::StrainTensor<float>			STPixelType;
  typedef itk::Image<STPixelType, Dimension>	STImageType;
  typedef itk::FixedArray<RealType,2>                   EigenValuesArrayType;
  typedef itk::Matrix<RealType,2,2>                     EigenVectorsMatrixType;
  typedef itk::FixedArray<RealType,2>                   EigenVectorType;
  typedef itk::FixedArray<float,2>			DeformPixelType;
  typedef itk::Image<DeformPixelType,4>			DeformImageType;

  typedef enum ColorModes{DEF,EIG0,EIG1,ST0,ST1,ST2} ColorModes;

  static vtkTensorGlyphStrain *New();
  vtkTensorGlyphStrain();
  ~vtkTensorGlyphStrain();


  void SetInput(STImageType::Pointer entrada)
    {this->input = entrada;};

  void SetDeformImage(DeformImageType::Pointer entrada)
    {this->deformImage = entrada;};

  vtkPolyData *GetOutput();
  vtkPolyData *GetDeformOutput();

  void SetInputPoints(vtkPoints *points)
    {this->inputPoints = points;};
  vtkPoints* GetInputPoints()
    {return this->inputPoints;};

  void SetGlyphType(int type)
    {this->GlyphType = type;};
  int GetGlyphType()
    {return this->GlyphType;};

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

  void SetPlanoZ(int value)
    {this->PlanoZ = value;};
  int GetPlanoZ()
    {return this->PlanoZ;};

  void SetTiempo(int value)
    {this->Tiempo = value;};
  int GetTiempo()
    {return this->Tiempo;};

  void SetDeformArray(vtkFloatArray *array)
    {this->deformArray = array;};
  vtkFloatArray* GetDeformArray()
    {return this->deformArray;};

  void interpolarTensor (double x[3], double eigval_out[2], double *angulo, double deform_values[2]);

protected:

  STImageType::Pointer input;
  DeformImageType::Pointer deformImage;
  vtkPoints *inputPoints;
  int GlyphType;

  int Scaling; // Determine whether scaling of geometry is performed
  double ScaleFactor; // Scale factor to use to scale geometry
  int ColorMode; // The coloring mode to use for the glyphs.

  int PlanoZ, Tiempo;

  vtkFloatArray *deformArray;
};

#endif
