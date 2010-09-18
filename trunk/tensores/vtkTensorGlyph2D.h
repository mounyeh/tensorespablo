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
  typedef itk::StrainTensor<float>			StrainPixelType;
  typedef itk::Image<StrainPixelType, Dimension>	StrainImageType;
  typedef itk::FixedArray<RealType,2>                   EigenValuesArrayType;
  typedef itk::Matrix<RealType,2,2>                     EigenVectorsMatrixType;
  typedef itk::FixedArray<RealType,2>                   EigenVectorType;

  typedef enum ColorModes{INV,EIG0,EIG1,ST0,ST1,ST2} ColorModes;

  static vtkTensorGlyphStrain *New();
  vtkTensorGlyphStrain();
  ~vtkTensorGlyphStrain();


  void SetInput(StrainImageType::Pointer entrada)
    {this->input = entrada;};

  vtkPolyData *GetOutput();
  vtkPolyData *GetDeformOutput();

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


  void ComputeEigenSystem(vtkIdType,double*[2],double[2]);

protected:

  StrainImageType::Pointer input;
  vtkPoints *inputPoints;
  int GlyphType;

  int Scaling; // Determine whether scaling of geometry is performed
  double ScaleFactor; // Scale factor to use to scale geometry
  int ColorMode; // The coloring mode to use for the glyphs.

  int PlanoZ, Tiempo;

};

#endif
