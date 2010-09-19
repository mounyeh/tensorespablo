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

  void SetScaling(bool value)
    {this->Scaling = value;};
  bool GetScaling()
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

  void interpolarTensor (double x[3], double eigval_out[2], double *angulo, double deform_values[2]);

protected:

  STImageType::Pointer input;		// Imagen tensorial de entrada
  DeformImageType::Pointer deformImage;	// Campo de deformaciones de entrada
  vtkPoints *inputPoints;		// Puntos de entrada, para mostrar glifos en puntos que no pertenecen a la imagen

  bool Scaling; 			// Determina si se aplica la escala
  double ScaleFactor;
  int ColorMode; 			// Tipo de coloreado para los glifos

  int PlanoZ;				// Corte a visualizar
  int Tiempo;				// Instante de tiempo a visualizar

};

#endif
