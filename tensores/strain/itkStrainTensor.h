#ifndef __itkStrainTensor_h
#define __itkStrainTensor_h

// Undefine an eventual StrainTensor macro
#ifdef StrainTensor
#undef StrainTensor
#endif

#include <itkFixedArray.h>
#include <itkMatrix.h>
#include <itkNumericTraits.h>
#include <itkRGBPixel.h>

//#define TOL 1e-10
#define TOL 1e-10

/*
Revisar: FA, RA, Mode, Deviatoric, mean diffusivity, shape coefficients, computelogvector, computefromlogvector, getrgbcolorcode
Hacer: eigenvalues, eigenvector, eigensystem, computetensor
*/
namespace itk
{

template < typename TComponent >
class StrainTensor: public FixedArray<TComponent,3>
{
public:
	/** Standard class typedefs. */
	typedef StrainTensor                Self;
	typedef FixedArray<TComponent,3> Superclass;

	/** Propagating some typedef from the superclass */
	typedef TComponent                                        ComponentType;
	typedef ComponentType                                     ComponentArrayType[6];
	typedef typename itk::NumericTraits<TComponent>::RealType RealType;
	typedef itk::FixedArray<RealType,2>                       EigenValuesArrayType;
	typedef itk::Matrix<RealType,2,2>                         EigenVectorsMatrixType;
	typedef itk::FixedArray<RealType,2>                       EigenVectorType;

	typedef itk::RGBPixel<TComponent>						  RGBPixelType;
	
	/** Default Constructor. */
	StrainTensor();

	/** Constructor with initialization. */
	StrainTensor(const Self& r);
	StrainTensor(const Superclass& r);
	StrainTensor(const ComponentType& r);
	StrainTensor(const ComponentArrayType r);

	/** Returns the number of components in this vector type */
	static unsigned int GetNumberOfComponents(){ return 6;}
	
	void SetNthComponent(int c, const ComponentType& v)  
    { 
		this->operator[](c) = v; 
	}

	/** Pass-through assignment operator for the Array base class. */
	Self& operator= (const Self& r);
	Self& operator= (const Superclass & r);
	Self& operator= (const ComponentType& r);
	Self& operator= (const ComponentArrayType r);
	
	/** Addition and substraction operators */
	Self operator+(const Self &tensor) const;
	Self operator-(const Self &tensor) const;
	const Self & operator+=(const Self &tensor);
	const Self & operator-=(const Self &tensor);
	
	/** Addition and substraction operators with scalars*/
	Self operator+(const RealType & scalar) const;
	Self operator-(const RealType & scalar) const;
	const Self & operator+=(const RealType & scalar);
	const Self & operator-=(const RealType & scalar);
	
	/** Product/division with scalars */
	Self operator*(const RealType & scalar ) const;
	Self operator/(const RealType & scalar ) const;
	const Self & operator*=(const RealType & scalar );
	const Self & operator/=(const RealType & scalar );
	
	/** Tensor comparisons */
	bool operator>( const Self &tensor ) const;
	bool operator>=( const Self &tensor ) const;
	bool operator<( const Self &tensor ) const;
	bool operator<=( const Self &tensor ) const;
	
	/** Scalar comparisons */
	bool operator>( const RealType& scalar ) const;
	bool operator>=( const RealType& scalar ) const;
	bool operator<( const RealType& scalar ) const;
	bool operator<=( const RealType& scalar ) const;
	
	/** Inner product between tensors */
	RealType operator^( const Self & tensor ) const;
	
	/** Matrix product between tensors */
	Self operator*( const Self & tensor ) const;
	const Self& operator*=( const Self & tensor );
	
	/** Matrix-based indexing */
	ComponentType & operator()( unsigned int row, unsigned int col );
	const ComponentType & operator()( unsigned int row, unsigned int col ) const;
	
	/** Identity matrix/tensor */
	void SetIdentity( void );
	
	/** Get Trace value */
	RealType GetTrace() const;

	/** Get the value of the determinant */
	RealType GetDeterminant() const;
	
	/** Get the mode of the tensor */
	RealType GetMode() const;
	
	/** Get the deviatoric of the tensor. */
	Self GetDeviatoric() const;

	/** Get the invariant (as described in Selskog et al.)	*/
	RealType GetInvariant() const;
	
	/** Compute the eigenvalues of the tensor */
	void ComputeEigenValues( EigenValuesArrayType& eig ) const;

	/** Compute one single eigenvector associated to a given eigen-vector */
	void ComputeEigenVector( RealType eigval, RealType* eigvec ) const;
	
	/** Compute the whole eigen-system */
	void ComputeEigenSystem( EigenValuesArrayType& eigval, EigenVectorsMatrixType& eigvec );
	
	/** Get the color code for RGB */ 
	RGBPixelType GetRGBColorCode();

	
	/** Compute the tensor from its eigen-system. Eigenvectors are arranged in columns */
	void ComputeTensor( EigenVectorsMatrixType eigvec, EigenValuesArrayType eigval );
	
	/** Compute the vector in the Log-Euclidean space */
	Self ComputeLogVector() const;
	/** Compute the vector in the Log-Euclidean space with weigthing factors*/
	Self ComputeLogVector(Vector<TComponent,6>) const;
	
	/** Compute the tensor from vector in the Log-Euclidean space */
	Self ComputeTensorFromLogVector() const;
	/** Compute the tensor from vector in the Log-Euclidean space with weighting factors  */
	Self ComputeTensorFromLogVector(Vector<TComponent,6>) const;

	
//public:
//	static const RealType m_Tolerance=1e-5;
};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_StrainTensor(_, EXPORT, x, y) namespace itk { \
  _(1(class EXPORT StrainTensor< ITK_TEMPLATE_1 x >)) \
  namespace Templates { typedef StrainTensor< ITK_TEMPLATE_1 x > \
                                            StrainTensor##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkStrainTensor+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkStrainTensor.txx"
#endif

#endif

