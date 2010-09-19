#ifndef _itkStrainTensor_txx
#define _itkStrainTensor_txx

#include "itkStrainTensor.h"
#include <itkSymmetricEigenAnalysis.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkMatrix.h>
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"

namespace itk
{
/*---------------------------------------------------------------------------**/
/** Constructors...*/
template<class T>
StrainTensor<T>::StrainTensor(){
	this->Fill( itk::NumericTraits<ComponentType>::Zero );
}

template<class T>
StrainTensor<T>::StrainTensor( const Self & r ):Superclass(r){
}

template<class T>
StrainTensor<T>::StrainTensor( const Superclass & r ):Superclass(r){
}


template<class T>
StrainTensor<T>::StrainTensor( const ComponentType & r )
{
	this->Fill( r );
}

template<class T>
StrainTensor<T>::StrainTensor( const ComponentArrayType r ):Superclass(r){
} 
/*---------------------------------------------------------------------------**/


/*---------------------------------------------------------------------------**/
/** Assignment operators...*/
template<class T>
StrainTensor<T>& StrainTensor<T>::operator= (const Self& r)
{
	Superclass::operator=(r);
	return *this;
}

template<class T>
StrainTensor<T>& StrainTensor<T>::operator= (const ComponentType & r)
{
	Superclass::operator=(r);
	return *this;
}

template<class T>
StrainTensor<T>& StrainTensor<T>::operator= (const ComponentArrayType r)
{
	Superclass::operator=(r);
	return *this;
}

template<class T>
StrainTensor<T>& StrainTensor<T>::operator= (const Superclass & r)
{
	Superclass::operator=(r);
	return *this;
}
/*---------------------------------------------------------------------------**/

/*---------------------------------------------------------------------------**/
/** Addition and substraction operators... */
template<class T>
StrainTensor<T> StrainTensor<T>::operator+(const StrainTensor<T> &tensor) const
{
	Self result;
	result[0] = (*this)[0] + tensor[0];
	result[1] = (*this)[1] + tensor[1];
	result[2] = (*this)[2] + tensor[2];
	return result;
}
template<class T>
StrainTensor<T> StrainTensor<T>::operator-(const StrainTensor<T> &tensor) const
{
	Self result;
	result[0] = (*this)[0] - tensor[0];
	result[1] = (*this)[1] - tensor[1];
	result[2] = (*this)[2] - tensor[2];
	return result;
}
template<class T>
const StrainTensor<T>& StrainTensor<T>::operator+=(const StrainTensor<T> &tensor)
{
	(*this)[0] += tensor[0];
	(*this)[1] += tensor[1];
	(*this)[2] += tensor[2];
	return *this;
}
template<class T>
const StrainTensor<T>& StrainTensor<T>::operator-=(const StrainTensor<T> &tensor)
{
	(*this)[0] -= tensor[0];
	(*this)[1] -= tensor[1];
	(*this)[2] -= tensor[2];
	return *this;
}
/*---------------------------------------------------------------------------**/

/*---------------------------------------------------------------------------**/
/** Addition and substraction operators with scalar values... */
template<class T>
StrainTensor<T> StrainTensor<T>::operator+(const RealType & scalar) const
{
	Self result;
	result[0] = (*this)[0] + scalar;
	result[1] = (*this)[1] + scalar;
	result[2] = (*this)[2] + scalar;
	return result;
}
template<class T>
StrainTensor<T> StrainTensor<T>::operator-(const RealType & scalar) const
{
	Self result;
	result[0] = (*this)[0] - scalar;
	result[1] = (*this)[1] - scalar;
	result[2] = (*this)[2] - scalar;
	return result;
}
template<class T>
const StrainTensor<T>& StrainTensor<T>::operator+=(const RealType & scalar)
{
	(*this)[0] += scalar;
	(*this)[1] += scalar;
	(*this)[2] += scalar;
	return *this;
}
template<class T>
const StrainTensor<T>& StrainTensor<T>::operator-=(const RealType & scalar)
{
	(*this)[0] -= scalar;
	(*this)[1] -= scalar;
	(*this)[2] -= scalar;
	return *this;
}
/*---------------------------------------------------------------------------**/

/*---------------------------------------------------------------------------**/
/** Product/division with scalars */
template<class T>
StrainTensor<T> StrainTensor<T>::operator*(const RealType &scalar) const
{
	Self result;
	result[0] = (*this)[0] * scalar;
	result[1] = (*this)[1] * scalar;
	result[2] = (*this)[2] * scalar;
	return result;
}
template<class T>
StrainTensor<T> StrainTensor<T>::operator/(const RealType &scalar) const
{
	Self result;
	RealType scalar2 = 1.0f/scalar;
	result[0] = (*this)[0] * scalar2;
	result[1] = (*this)[1] * scalar2;
	result[2] = (*this)[2] * scalar2;
	return result;
}
template<class T>
const StrainTensor<T>& StrainTensor<T>::operator*=(const RealType &scalar)
{
	(*this)[0] *= scalar;
	(*this)[1] *= scalar;
	(*this)[2] *= scalar;
	return *this;
}
template<class T>
const StrainTensor<T>& StrainTensor<T>::operator/=(const RealType &scalar)
{
	RealType scalar2 = 1.0f/scalar;
	(*this)[0] *= scalar2;
	(*this)[1] *= scalar2;
	(*this)[2] *= scalar2;
	return *this;
}
/*---------------------------------------------------------------------------**/
	
/*---------------------------------------------------------------------------**/
/** Matrix-based indexing. No bounds checking!!! */
template<class T>
const typename StrainTensor<T>::ComponentType&
StrainTensor<T>::operator()( unsigned int row, unsigned int col ) const
{
	return ((*this)[row+col]);
}
template<class T>
typename StrainTensor<T>::ComponentType&
StrainTensor<T>::operator()( unsigned int row, unsigned int col )
{
	return ((*this)[row+col]);
}
/*---------------------------------------------------------------------------**/

	
/*---------------------------------------------------------------------------**/
/** Tensor comparisons */
template<class T>
bool StrainTensor<T>::operator>( const Self& tensor ) const
{
	return (  this->GetTrace() > tensor.GetTrace()  );
}
template<class T>
bool StrainTensor<T>::operator>=( const Self& tensor ) const
{
	return (  this->GetTrace() >= tensor.GetTrace()  );
}
template<class T>
bool StrainTensor<T>::operator<( const Self& tensor ) const
{
	return (  this->GetTrace() < tensor.GetTrace()  );
}
template<class T>
bool StrainTensor<T>::operator<=( const Self& tensor ) const
{
	return (  this->GetTrace() <= tensor.GetTrace()  );
}
/*---------------------------------------------------------------------------**/
	

	
/*---------------------------------------------------------------------------**/
/** Scalar comparisons */
template<class T>
bool StrainTensor<T>::operator>( const RealType& scalar ) const
{
	return (  this->GetTrace() > scalar  );
}
template<class T>
bool StrainTensor<T>::operator>=( const RealType& scalar ) const
{
	return (  this->GetTrace() >= scalar  );
}
template<class T>
bool StrainTensor<T>::operator<( const RealType& scalar ) const
{
	return (  this->GetTrace() < scalar  );
}
template<class T>
bool StrainTensor<T>::operator<=( const RealType& scalar ) const
{
	return (  this->GetTrace() <= scalar  );
}
/*---------------------------------------------------------------------------**/
	
	
	
/** Inner product between tensors */
template<class T>
typename StrainTensor<T>::RealType
StrainTensor<T>::operator^( const StrainTensor<T> & tensor ) const
{
	RealType product = itk::NumericTraits<RealType>::Zero;
	product += (*this)[0] * tensor[0];
	product += 2 * (*this)[1] * tensor[1];
	product += (*this)[2] * tensor[2];
	return product;
}


/*---------------------------------------------------------------------------**/
/** Matrix product between tensors */
template<class T>
StrainTensor<T> StrainTensor<T>::operator*( const StrainTensor<T> & tensor ) const
{
	Self result;
	result[0] = ((*this)[0])*tensor[0] + ((*this)[1])*tensor[1];
	result[1] = ((*this)[0])*tensor[1] + ((*this)[1])*tensor[2];
	result[2] = ((*this)[1])*tensor[1] + ((*this)[2])*tensor[2];
	return result;
}

template<class T>
const StrainTensor<T>& StrainTensor<T>::operator*=( const StrainTensor<T> & tensor )
{
	Self result;
	result[0] = ((*this)[0])*tensor[0] + ((*this)[1])*tensor[1];
	result[1] = ((*this)[0])*tensor[1] + ((*this)[1])*tensor[2];
	result[2] = ((*this)[1])*tensor[1] + ((*this)[2])*tensor[2];
	(*this) = result;
	return (*this);
}
/*---------------------------------------------------------------------------**/
	
	
/** Identity matrix/tensor */
template<class T>
void StrainTensor<T>::SetIdentity()
{
	this->Fill( itk::NumericTraits<ComponentType>::Zero );
	(*this)[0] = static_cast<ComponentType>( 1 );
	(*this)[2] = static_cast<ComponentType>( 1 );
}

/** Get the trace */
template<class T>
typename StrainTensor<T>::RealType
StrainTensor<T>::GetTrace() const
{
	return (   static_cast<RealType>( (*this)[0] )   +   static_cast<RealType>( (*this)[2] )   );
}

	
/**
 *  Compute the value of relative anisotropy
 */
template<class T>
typename StrainTensor<T>::RealType
StrainTensor<T>::GetDeterminant() const
{
	RealType det = itk::NumericTraits<RealType>::Zero;
	
	det += ((*this)[0])*((*this)[2]);
	det -= ((*this)[1])*((*this)[1]);
	
	return det;
}


/**
 *  Compute the invariant (the Frobenius norm, as described in Selskog et al.)
 */
template<class T>
typename StrainTensor<T>::RealType
StrainTensor<T>::GetInvariant() const
{
	EigenValuesArrayType eig;
	this->ComputeEigenValues( eig );
	RealType invariant = eig[0]*eig[0] + eig[1]*eig[1];
	return ( invariant );
}
		
	
/** Compute the eigenvalues of the tensor. Direct method!!! */
template<class T>
void StrainTensor<T>::ComputeEigenValues( EigenValuesArrayType& eig ) const
{

	if(fabs((*this)[1]*(*this)[1]) <= 0.1e-20 ) {
		eig[0] = (*this)[0];
		eig[1] = (*this)[2];
		return;
	}

	// First of all, compute the eigenvalues, together with the rank of the filter:
	double tr = (*this)[0] + (*this)[2];
	double det = (*this)[0] * (*this)[2] - (*this)[1] * (*this)[1];
	double S = sqrt( (tr/2)*(tr/2) - det );
	eig[0] = tr/2 + S;
	eig[1] = tr/2 - S;


	return;
}

	

/** Compute the whole eigen-system */
template<class T>
void StrainTensor<T>
::ComputeEigenSystem( EigenValuesArrayType& eigval, EigenVectorsMatrixType& eigvec )
{

	if(fabs((*this)[1]*(*this)[1]) <= 0.1e-20 ) {
		eigval[0] = (*this)[0]; eigvec[0][0] = 1; eigvec[1][0] = 0;
		eigval[1] = (*this)[2]; eigvec[0][1] = 0; eigvec[1][1] = 1;
		return;
	}

	// First of all, compute the eigenvalues, together with the rank of the filter:

	double tr = (*this)[0] + (*this)[2];
	double det = (*this)[0] * (*this)[2] - (*this)[1] * (*this)[1];
	double S = sqrt( (tr/2)*(tr/2) - det );
	eigval[0] = tr/2 + S;
	eigval[1] = tr/2 - S;

	if ((*this)[1]!=0) {
		eigvec[0][0]=eigval[0]-(*this)[2];
		eigvec[1][0]=(*this)[1];
		eigvec[0][1]=eigval[1]-(*this)[2];
		eigvec[1][1]=(*this)[1];
	}

	double n1 = sqrt(eigvec[0][0]*eigvec[0][0]+eigvec[1][0]*eigvec[1][0]);
	eigvec[0][0] /= n1; eigvec[1][0] /= n1;

	double n2 = sqrt(eigvec[0][1]*eigvec[0][1]+eigvec[1][1]*eigvec[1][1]);
	eigvec[0][1] /= n2; eigvec[1][1] /= n2;

	return;

}
	

/** Compute the tensor from its eigen-system */
template<class T>
void StrainTensor<T>::ComputeTensor( EigenVectorsMatrixType eigvec, EigenValuesArrayType eigval )
{

	(*this)[0] = eigval[0]*eigvec[0][0]*eigvec[0][0] + eigval[1]*eigvec[0][1]*eigvec[0][1];
	(*this)[0] = eigval[0]*eigvec[0][0]*eigvec[1][0] + eigval[1]*eigvec[0][1]*eigvec[1][1];
	(*this)[0] = eigval[0]*eigvec[1][0]*eigvec[1][0] + eigval[1]*eigvec[1][1]*eigvec[1][1];

}

} // end namespace itk

#endif

