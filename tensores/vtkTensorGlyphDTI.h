#ifndef __vtkTensorGlyphDTI_h
#define __vtkTensorGlyphDTI_h

#include "tensor/itkDTITensor.h"
#include <itkFixedArray.h>
#include <itkMatrix.h>
#include <itkNumericTraits.h>

#include "tensor/itkDTITensor.h"
#include "tensor/itkDWImages.h"

#include "vtkPolyData.h"

class vtkTensorGlyphDTI
{

public:

	enum { Dimension = 3 };

	typedef itk::NumericTraits<float>::RealType 		RealType;
	typedef itk::DTITensor<float>				TensorPixelType;
	typedef itk::Image<TensorPixelType, Dimension>		TensorImageType;
	typedef itk::FixedArray<RealType,3>	 		EigenValuesArrayType;
	typedef itk::Matrix<RealType,3,3>	 		EigenVectorsMatrixType;
	typedef itk::FixedArray<RealType,3>	 		EigenVectorType;


	static vtkTensorGlyphDTI *New();
	vtkTensorGlyphDTI();
	~vtkTensorGlyphDTI();


	void GetOutput(vtkPolyData*);

	void interpolacionLineal (double x[3],TensorPixelType *pixel);
	void interpolacionLogEuclidea (double x[3],TensorPixelType *pixel);


	void SetInput(TensorImageType::Pointer image)
		{this->input = image;};
	TensorImageType::Pointer SetInput()
		{return this->input;};

	void SetInputPoints(vtkPoints *points)
		{this->inputPoints = points;};
	vtkPoints *GetInputPoints()
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

	void SetPhiResolution(int value)
		{this->PhiResolution = value;};
	int GetPhiResolution()
		{return this->PhiResolution;};

	void SetThetaResolution(int value)
		{this->ThetaResolution = value;};
	int GetThetaResolution()
		{return this->ThetaResolution;};

	void SetGamma(double value)
		{this->Gamma = value;};
	double GetGamma()
		{return this->Gamma;};

	void SetFilterMode(int value)
		{this->FilterMode = value;};
	int GetFilterMode()
		{return this->FilterMode;};

	void SetFilterThreshold(double value)
		{this->FilterThreshold = value;};
	double GetcsThreshold()
		{return this->FilterThreshold;};

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

	enum
	{
	ELLIPSOID,
	CUBOID,
	SUPERQUADRIC
	};

	enum
	{
	COLOR_BY_FA,
	COLOR_BY_RA,
	COLOR_BY_CL
	};

	enum
	{
	FILTER_BY_FA,
	FILTER_BY_CL,
	FILTER_BY_CS
	};

protected:

	TensorImageType::Pointer input;	// Imagen tensorial de entrada
	vtkPolyData *source;
	vtkPoints *inputPoints;		// Puntos de entrada donde mostrar glifos
	int GlyphType;

	int Scaling; 			// Determina si se aplica la escala
	double ScaleFactor; 
	int ColorMode; 			// Tipo de coloreado para los glifos

	int Bounds[6];			// Límites de la imagen en los que representar glifos

	int PhiResolution;
	int ThetaResolution;

	double Gamma;			// Parámetro gamma para supercuádricas

	int FilterMode;			// Magnitud a considerar para la discriminación de glifos
	double FilterThreshold;		// Umbral para la discriminación

};

#endif
