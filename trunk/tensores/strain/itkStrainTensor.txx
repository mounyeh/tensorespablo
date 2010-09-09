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
	//m_Tolerance = 1e-5;
}

template<class T>
StrainTensor<T>::StrainTensor( const Self & r ):Superclass(r){//m_Tolerance = 1e-5;
}

template<class T>
StrainTensor<T>::StrainTensor( const Superclass & r ):Superclass(r){//m_Tolerance = 1e-5;
}


template<class T>
StrainTensor<T>::StrainTensor( const ComponentType & r )
{
	this->Fill( r );
//	m_Tolerance = 1e-5;
}

template<class T>
StrainTensor<T>::StrainTensor( const ComponentArrayType r ):Superclass(r){//m_Tolerance = 1e-5;
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
 *  Compute the value of fractional anisotropy
 */
template<class T>
typename StrainTensor<T>::RealType
StrainTensor<T>::GetFractionalAnisotropy() const
{
	// Compute the squared Frobenius norm:
	RealType isp   = (*this)^(*this);
	if( isp > TOL ){
		RealType trace      = this->GetTrace();
		RealType anisotropy = ( 3*isp - trace*trace )/( 2*isp );
		return static_cast< RealType >( vcl_sqrt(anisotropy) );
    }
    return itk::NumericTraits<RealType>::Zero;
}


/**
 *  Compute the value of relative anisotropy
 */
template<class T>
typename StrainTensor<T>::RealType
StrainTensor<T>::GetRelativeAnisotropy() const
{
	// Basser and Perpaoli implementation, eq. [14]
	// Compute the trace and the Frobenius norm:
	RealType trace = this->GetTrace();
	RealType isp   = (*this)^(*this);
	
	// Avoid negative trace and traces small enough to look like a division by zero.
	if( trace < TOL )
		return itk::NumericTraits<RealType>::Zero;
	
	RealType anisotropy = 3*isp - trace*trace;
	
	if( anisotropy  < TOL )
		return itk::NumericTraits<RealType>::Zero;
	
	anisotropy = static_cast<RealType>(  vcl_sqrt(anisotropy )   );
	return (anisotropy/trace);
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
 *  Compute the mode of the tensor as in Kindlmann et al.
 */
template<class T>
typename StrainTensor<T>::RealType
StrainTensor<T>::GetMode() const
{
	RealType trace = this->GetTrace();
	RealType dNorm = (*this)^(*this) - (1.0f/3.0f)*trace*trace;
	if( dNorm < TOL )
		return itk::NumericTraits<RealType>::Zero;
	dNorm = sqrt( dNorm );
	dNorm = 7.3484692283495342946 / ( dNorm*dNorm*dNorm );
	return ( this->GetDeviatoric->GetDeterminant() * dNorm );
}
	
	
/**
 *  Compute the deviatoric tensor
 */
template<class T>
StrainTensor<T> StrainTensor<T>::GetDeviatoric() const
{
	Self identity;
	identity->SetIdentity();
	identity *= ( (1.0f/3.0f)*(this->GetTrace()) );
	return ( *this - identity );
}
	
/**
 *  Compute the mean diffusivity
*/
template<class T>
typename StrainTensor<T>::RealType
StrainTensor<T>::GetMeanDiffusivity() const
{
	return ( (this->GetTrace()/3));
}

/**
 *  Compute the invariant (as described in Selskog et al.)
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
		

/** Compute the geometric parameters of the tensor that describes its shape */
template<class T>
void StrainTensor<T>::ComputeShapeCoefficients( RealType& cl, RealType& cp, RealType& cs )
{
	cl = cp = cs = itk::NumericTraits<RealType>::Zero;
	// Compute the eigen-values:
	EigenValuesArrayType eig;
	this->ComputeEigenValues( eig );
	// Make sure the first eigenvalue is not null:
	if( eig[0] < TOL )
		return;
	// Inverse of the first eigen-value:
	RealType div = itk::NumericTraits<RealType>::One/eig[0];
	// Computation of shape parameters:
	cs  = eig[2] * div;
	cp  = eig[1] * div;
	cl  = 1 - cp;
	cp -= cs;
	return;
}
	
	
	
/** Compute the eigenvalues of the tensor. Direct method!!! */
template<class T>
void StrainTensor<T>::ComputeEigenValues( EigenValuesArrayType& eig ) const
{

    	if((*this)[1]*(*this)[1] <= 0.1e-20 ) {
    	    eig[0] = (*this)[0];
    	    eig[1] = (*this)[2];
    	    return;
    	}

    	RealType tr = (*this)[0] + (*this)[2];
    	RealType det = (*this)[0] * (*this)[2] - (*this)[1] * (*this)[1];
    	RealType S = sqrt( (tr/2)*(tr/2) - det );
    	eig[0] = tr/2 + S;
    	eig[1] = tr/2 - S;
/*
    	double SS = sqrt( std::max(square(((*this)[0]-(*this)[2])/2) + (*this)[1] * (*this)[1], 0.0) );
    	if( (*this)[0] - (*this)[2] < 0 ) {
        	*v1x = (*this)[1];
        	*v1y = - ((*this)[0]-(*this)[2])/2 + SS;
        	*v2x = + ((*this)[0]-(*this)[2])/2 - SS;
        	*v2y = (*this)[1];
    	} else {
        	*v2x = (*this)[1];
        	*v2y = - ((*this)[0]-(*this)[2])/2 - SS;
        	*v1x = + ((*this)[0]-(*this)[2])/2 + SS;
        	*v1y = (*this)[1];
    	}

    	double n1 = sqrt(square(*v1x)+square(*v1y));
    	*v1x /= n1; *v1y /= n1;
    	double n2 = sqrt(square(*v2x)+square(*v2y));
	*v2x /= n2; *v2y /= n2;
*/

/*	// Definition of the coefficients of the characteristic 
	// polynomial (the coefficient in lambda^3 is simply 1):
	
	RealType d12s = ((*this)[1])*((*this)[1]);
	RealType d13s = ((*this)[2])*((*this)[2]);
	RealType d23s = ((*this)[4])*((*this)[4]);
	
	RealType b  =  -(   ((*this)[0]) + ((*this)[3]) + ((*this)[5])   );
	RealType c  =  ((*this)[0])*((*this)[3]) + ((*this)[0])*((*this)[5]) + ((*this)[3])*((*this)[5]);
	c          -=  ( d12s + d13s + d23s );
	RealType d  =  -((*this)[0])*((*this)[3])*((*this)[5]) - 2*((*this)[1])*((*this)[2])*((*this)[4]);
	d          +=  d13s*((*this)[3]) + d12s*((*this)[5]) + d23s*((*this)[0]);
	
	// Computation of auxiliar values:
	RealType disc = 12*c*c*c - 3*b*b*c*c - 54*b*c*d + 81*d*d + 12*b*b*b*d;
	RealType alpha = c - (1.0f/3.0f)*b*b;
	RealType res  = (1.0f/3.0f)*b;
	RealType xiR = 36*b*c - 108*d - 8*b*b*b;
	RealType xiI = itk::NumericTraits<RealType>::Zero;
	RealType lambda1;
	RealType lambda2;
	RealType lambda3;
	if( disc >= itk::NumericTraits<RealType>::Zero ){
		RealType xiRi;
		xiR     += 12*::sqrt(  disc );
		if( xiR > TOL ){
			xiR      = ::pow( xiR, (double)(1.0f/3.0f) );
			xiRi = 1.0f/xiR;
		}
		else if( xiR < -TOL ){
			xiR      = -::pow( -xiR, (double)(1.0f/3.0f) );
			xiRi = 1.0f/xiR;
		}
		else
			xiR = xiRi = itk::NumericTraits<RealType>::Zero;
		lambda1  =  (1.0f/6.0f)*xiR - 2*alpha*xiRi;
		lambda2  = -0.5f*lambda1 - res;
		lambda1 -= res;
		lambda3  = lambda2;
	}
	else{
		xiI  = 12*::sqrt( -disc );
		RealType phi = static_cast<RealType>( 0.523598775598299 );
		if(   vcl_abs(xiR) > TOL   )
			phi = (1.0f/3.0f)*::atan( xiI/xiR );
		RealType xiAbs = ::pow( xiI*xiI+xiR*xiR, (double)(1.0f/6.0f) );
		RealType cphi  = ::cos( phi );
		RealType sphi  = ::sin( phi );
		xiR = xiAbs*cphi;
		xiI = xiAbs*sphi;
		RealType xiAbsi = 1.0f/xiAbs;
		RealType xiRi   =  xiAbsi*cphi;
		RealType xiIi   = -xiAbsi*sphi;
		
		lambda1 = (1.0f/6.0f)*xiR - 2*alpha*xiRi;
		lambda2 = -0.5f*lambda1;
		lambda3 = lambda2;
		RealType imPart = (0.144337567297406*xiI + 1.732050807568877*xiIi*alpha);
		lambda2 -= imPart;
		lambda3 += imPart;
		lambda1 -= res;
		lambda2 -= res;
		lambda3 -= res;
	}
	// Order the eigenvalues:
	eig[0]=lambda1; eig[1]=lambda2; eig[2]=lambda3;
	if( lambda1<lambda3 ){
		eig[0] = lambda3;
		eig[2] = lambda1;
	}
	if( lambda2<eig[2] ){
		eig[1] = eig[2];
		eig[2] = lambda2;
	}
	else if( lambda2>eig[0] ){
		eig[1] = eig[0];
		eig[0] = lambda2;
	}
*/
	return;
}

	
/** Compute one single eigenvector associated to a given eigen-vector */
// Be careful!!! This method may only be used once you have checked that
// eigval is an eigenvalue with multiplicity 1!!!!
template<class T>
void StrainTensor<T>
::ComputeEigenVector( RealType eigval, RealType* eigvec ) const
{
	// Compute an adequate tolerance:
	double e_tol = ( (*this)[0] + (*this)[3] + (*this)[5] );
	e_tol        = 0.001f * e_tol * e_tol;
	// Compute the corrected matrix:
	RealType matrix[3][3];
	matrix[0][0] = (*this)[0] - eigval;
	matrix[1][1] = (*this)[3] - eigval;
	matrix[2][2] = (*this)[5] - eigval;
	matrix[0][1] = matrix[1][0] = (*this)[1];
	matrix[0][2] = matrix[2][0] = (*this)[2];
	matrix[1][2] = matrix[2][1] = (*this)[4];
	// Find a minor with non-null determinant; there must be at least one, since
	// the eigen-value has multiplicity 1
	for( unsigned int r1=0; r1<2; ++r1 ){ // The whole four nested loops comprise at most 9 computations of det
		for( unsigned int r2=r1+1; r2<3; ++r2 ){
			for( unsigned int c1=0; c1<2; ++c1 ){
				for( unsigned int c2=c1+1; c2<3; ++c2 ){
					RealType det = matrix[r1][c1]*matrix[r2][c2] - matrix[r2][c1]*matrix[r1][c2];
					if( vcl_abs(det) > e_tol ){
						// We have found the minor, find the other column:
						unsigned int cr = 1 + (c1+c2)%2 - c1 - c1;
						// Find the independent term
						RealType b1 = -matrix[r1][cr];
						RealType b2 = -matrix[r2][cr];
						// Solve the system:
						det = 1/det;
						eigvec[cr] = 1;
						eigvec[c1] = det*( matrix[r2][c2]*b1 - matrix[r1][c2]*b2 );
						eigvec[c2] = det*( matrix[r1][c1]*b2 - matrix[r2][c1]*b1 );
						// Normalise to have unit norm:
						RealType norm = eigvec[0]*eigvec[0] + eigvec[1]*eigvec[1] + eigvec[2]*eigvec[2];
						norm = 1/vcl_sqrt( norm );
						eigvec[0] *= norm;
						eigvec[1] *= norm;
						eigvec[2] *= norm;
						// Return!
						return;
					}
				}
			}
		}
	}
}

	
/** Compute the whole eigen-system */
template<class T>
void StrainTensor<T>
::ComputeEigenSystem( EigenValuesArrayType& eigval, EigenVectorsMatrixType& eigvec )
{

  if((*this)[1]*(*this)[2] <= 0.1e-20 ) {
    eigval[0] = (*this)[0]; eigvec[0][0] = 1; eigvec[1][0] = 0;
    eigval[1] = (*this)[3]; eigvec[0][1] = 0; eigvec[1][1] = 1;
    return;
  }

  // First of all, compute the eigenvalues, together with the rank of the filter:
  double tr = (*this)[0] + (*this)[3];
  double det = (*this)[0] * (*this)[3] - (*this)[1] * (*this)[2];
  double S = sqrt( (tr/2)*(tr/2) - det );
  eigval[0] = tr/2 + S;
  eigval[1] = tr/2 - S;
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

  if ((*this)[1]!=0) {
    eigvec[0][0]=eigval[0]-(*this)[3];
    eigvec[1][0]=(*this)[1];
    eigvec[0][1]=eigval[1]-(*this)[3];
    eigvec[1][1]=(*this)[1];
  }

  double n1 = sqrt(eigvec[0][0]*eigvec[0][0]+eigvec[1][0]*eigvec[1][0]);
  eigvec[0][0] /= n1; eigvec[1][0] /= n1;

  double n2 = sqrt(eigvec[0][1]*eigvec[0][1]+eigvec[1][1]*eigvec[1][1]);
  eigvec[0][1] /= n2; eigvec[1][1] /= n2;

  return;


/*
    	if((*this)[1]*(*this)[1] <= 0.1e-20 ) {
    	    eigval[0] = (*this)[0]; eigvec[0][0] = 1; eigvec[1][0] = 0;
    	    eigval[1] = (*this)[2]; eigvec[0][1] = 0; eigvec[1][1] = 1;
    	    return;
    	}

	// First of all, compute the eigenvalues, together with the rank of the filter:
	this->ComputeEigenValues( eigval );

    	double SS = sqrt( std::max(square(((*this)[0]-(*this)[2])/2) + (*this)[1] * (*this)[1], 0.0) );
    	if( (*this)[0] - (*this)[2] < 0 ) {
        	eigvec[0][0] = (*this)[1];
        	eigvec[1][0] = - ((*this)[0]-(*this)[2])/2 + SS;
        	eigvec[0][1] = + ((*this)[0]-(*this)[2])/2 - SS;
        	eigvec[1][1] = (*this)[1];
    	} else {
        	eigvec[0][1] = (*this)[1];
        	eigvec[1][1] = - ((*this)[0]-(*this)[2])/2 - SS;
        	eigvec[0][0] = + ((*this)[0]-(*this)[2])/2 + SS;
        	eigvec[1][0] = (*this)[1];
    	}

    	double n1 = sqrt((eigvec[0][0])*(eigvec[0][0])+(eigvec[1][0])*(eigvec[1][0]));
    	eigvec[0][0] /= n1; eigvec[1][0] /= n1;
    	double n2 = sqrt((eigvec[0][1])*(eigvec[0][1])+(eigvec[1][1])*(eigvec[1][1]));
	eigvec[0][1] /= n2; eigvec[1][1] /= n2;

	return;
*/
}
	
/** Get color code: RGB value coded the major eigenvector weighted by FA **/
template<class T>                             
typename StrainTensor<T>::RGBPixelType           
StrainTensor<T>::GetRGBColorCode(){              
    EigenValuesArrayType eigval;              
	ComputeEigenValues(eigval);                  
	                                             
	RealType eigvec[3];                          
	ComputeEigenVector( eigval[0], eigvec);      
                                              
    RealType fa;                              
	fa=GetFractionalAnisotropy();                
	                                             
	RGBPixelType colorCoding;                    
                                              
	for(unsigned i=0; i<3; i++){                 
		colorCoding[i]=fa*fabs(eigvec[i]);          
	}                                            
	return colorCoding;
}

/** Esta es otra forma, calculando los autovectores y autovalores como lo hace itk **/
/** Get color code: RGB value coded the major eignevector weighted by FA **/
/*
template<class T>
typename StrainTensor<T>::RGBPixelType
StrainTensor<T>::GetRGBColorCode(){    

    itk::SymmetricSecondRankTensor<float,3> symTensor;
	for(unsigned i=0; i<6; i++){
		symTensor[i]=(*this)[i];
	}
	
	itk::SymmetricSecondRankTensor<float,3>::EigenValuesArrayType eigenValues;
	itk::SymmetricSecondRankTensor<float,3>::EigenVectorsMatrixType eigenVectors;
	//ComputeEigenValues( eigenValues );
	symTensor.ComputeEigenAnalysis(eigenValues, eigenVectors);
    
	RealType fa;
	fa=GetFractionalAnisotropy();
    RGBPixelType colorCoding;   
	for(unsigned i=0; i<3; i++){
		colorCoding[i]=fa*fabs(eigenVectors[2][i]);
	}
	return colorCoding;
}
*/
	
/** Compute the tensor from its eigen-system */
template<class T>
void StrainTensor<T>::ComputeTensor( EigenVectorsMatrixType eigvec, EigenValuesArrayType eigval )
{
	(*this)[0] = eigval[0]*eigvec(0,0)*eigvec(0,0) + eigval[1]*eigvec(0,1)*eigvec(0,1) + eigval[2]*eigvec(0,2)*eigvec(0,2);
	(*this)[1] = eigval[0]*eigvec(0,0)*eigvec(1,0) + eigval[1]*eigvec(0,1)*eigvec(1,1) + eigval[2]*eigvec(0,2)*eigvec(1,2);
	(*this)[2] = eigval[0]*eigvec(0,0)*eigvec(2,0) + eigval[1]*eigvec(0,1)*eigvec(2,1) + eigval[2]*eigvec(0,2)*eigvec(2,2);
	(*this)[3] = eigval[0]*eigvec(1,0)*eigvec(1,0) + eigval[1]*eigvec(1,1)*eigvec(1,1) + eigval[2]*eigvec(1,2)*eigvec(1,2);
	(*this)[4] = eigval[0]*eigvec(1,0)*eigvec(2,0) + eigval[1]*eigvec(1,1)*eigvec(2,1) + eigval[2]*eigvec(1,2)*eigvec(2,2);
	(*this)[5] = eigval[0]*eigvec(2,0)*eigvec(2,0) + eigval[1]*eigvec(2,1)*eigvec(2,1) + eigval[2]*eigvec(2,2)*eigvec(2,2);
}


/**
* Compute the vector in the Log-Euclidean space
*/
template<class T>
StrainTensor<T>
StrainTensor<T>
::ComputeLogVector() const
{
	EigenValuesArrayType eigval;
	EigenVectorsMatrixType eigvec;
	ComputeEigenSystem(eigval, eigvec);
	
	Matrix<T, 3, 3> logEigval;
	logEigval.Fill(0.0);
	
	for(unsigned i=0; i<3; i++){
		logEigval(i,i)=log(eigval(i));
			} 	
	Matrix<T,3,3> logMatrix;
	logMatrix=eigvec.GetTranspose()*logEigval*eigvec;
	
	StrainTensor<T> result;
	result[0]=logMatrix[0][0];
	result[1]=logMatrix[1][1];	
	result[2]=logMatrix[2][2];
	result[3]=sqrt(2)*logMatrix[0][1];
	result[4]=sqrt(2)*logMatrix[0][2];
	result[5]=sqrt(2)*logMatrix[1][2];
	
	return result;
}


/**
* Compute the vector in the Log-Euclidean space with weigthing factor
*/
template<class T>
StrainTensor<T>
StrainTensor<T>
::ComputeLogVector(Vector<T,6> weights) const
{
	EigenValuesArrayType eigval;
	EigenVectorsMatrixType eigvec;
	ComputeEigenSystem(eigval, eigvec);
	
	Matrix<T, 3, 3> logEigenValues;
	logEigenValues.Fill(0.0);
	
	for(unsigned i=0; i<3; i++){
		logEigenValues(i,i)=log(eigval(i));
			} 	
	Matrix<T,3,3> logMatrix;
	logMatrix=eigvec.GetTranspose()*logEigenValues*eigvec;
	
	StrainTensor<T> result;
	result[0]=weights[0]*logMatrix[0][0];
	result[1]=weights[1]*logMatrix[1][1];	
	result[2]=weights[2]*logMatrix[2][2];
	result[3]=weights[3]*logMatrix[0][1];
	result[4]=weights[4]*logMatrix[0][2];
	result[5]=weights[5]*logMatrix[1][2];
	
	return result;
}

/**
* Translate the tensor to the original space from the log euclidean space
*/
template<class T>
StrainTensor<T>
StrainTensor<T>
::ComputeTensorFromLogVector() const
{
	EigenValuesArrayType logEigenValues;
	EigenVectorsMatrixType logEigenVectors;

	this[3]/=sqrt(2);
	this[4]/=sqrt(2);
	this[5]/=sqrt(2);

	ComputeEigenSystem(logEigenValues, logEigenVectors);
	
	Matrix<T, 3, 3> eigval;
	eigval.Fill(0.0);
	
	for(unsigned i=0; i<3; i++){
		eigval(i,i)=exp(logEigenValues(i));
			} 	
			
	Matrix<T,3,3> tensorMatrix;
	tensorMatrix=logEigenVectors.GetTranspose()*eigval*logEigenVectors;
	
	StrainTensor<T> result;	
	result[0]=tensorMatrix[0][0];
	result[1]=tensorMatrix[0][1];	
	result[2]=tensorMatrix[0][2];
	result[3]=tensorMatrix[1][1];
	result[4]=tensorMatrix[1][2];
	result[5]=tensorMatrix[2][2];
	
	return result;
}

/**
* Translate the tensor to the original space from the log euclidean space
*/
template<class T>
StrainTensor<T>
StrainTensor<T>
::ComputeTensorFromLogVector(Vector<T,6> weights) const
{
	EigenValuesArrayType logEigenValues;
	EigenVectorsMatrixType logEigenVectors;

	for(unsigned i=0; i<6; i++){
		this[i]/=weights[i];
	}
	ComputeEigenSystem(logEigenValues, logEigenVectors);
	
	Matrix<T, 3, 3> eigval;
	eigval.Fill(0.0);
	
	for(unsigned i=0; i<3; i++){
		eigval(i,i)=exp(logEigenValues(i));
			} 	
			
	Matrix<T,3,3> tensorMatrix;
	tensorMatrix=logEigenVectors.GetTranspose()*eigval*logEigenVectors;
	
	StrainTensor<T> result;	
	result[0]=tensorMatrix[0][0];
	result[1]=tensorMatrix[0][1];	
	result[2]=tensorMatrix[0][2];
	result[3]=tensorMatrix[1][1];
	result[4]=tensorMatrix[1][2];
	result[5]=tensorMatrix[2][2];
	
	return result;
}



} // end namespace itk

#endif

