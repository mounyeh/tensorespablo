/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkD3BSplineTransform.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:29:28 $
  Version:   $Revision: 1.29 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#ifndef __itkD3BSplineTransform_h
#define __itkD3BSplineTransform_h

//#include <iostream>
#include "itkTransform.h"
#include "itkImage.h"
#include "itkExceptionObject.h"
//#include "itkMatrix.h"

namespace itk
{

/** \brief D3BSpline transformation of a vector space (e.g. space coordinates)
 *
 * The same functionality could be obtained by using the Affine tranform,
 * but with a large difference in performace.
 *
 * \ingroup Transforms
 */
template < class TScalarType=double,          // Data type for scalars (float or double)
           unsigned int NDimensions=3>        // Number of dimensions
class ITK_EXPORT D3BSplineTransform : public Transform< TScalarType, NDimensions, NDimensions >
{
public:
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  /** Standard class typedefs. */
  typedef D3BSplineTransform Self;
  typedef Transform< TScalarType, NDimensions, NDimensions > Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------



  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  /** Image-related typedefs. Note that origin and spacing are of type double. */
  typedef typename Image<TScalarType,NDimensions>::IndexType IndexType;
  typedef typename Image<TScalarType,NDimensions>::SpacingType SpacingType;
  typedef typename Image<TScalarType,NDimensions>::SizeType SizeType;
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------


  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( D3BSplineTransform, Transform );

  /** Dimension of the domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------


  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard parameters container. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Standard Jacobian container. */
  typedef typename Superclass::JacobianType JacobianType;
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------


  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  /** Standard vector type for this class. */
  typedef Vector<TScalarType, itkGetStaticConstMacro(SpaceDimension)> InputVectorType;
  typedef Vector<TScalarType, itkGetStaticConstMacro(SpaceDimension)> OutputVectorType;

  /** Standard covariant vector type for this class. */
  typedef CovariantVector<TScalarType, itkGetStaticConstMacro(SpaceDimension)> InputCovariantVectorType;
  typedef CovariantVector<TScalarType, itkGetStaticConstMacro(SpaceDimension)> OutputCovariantVectorType;
  
  /** Standard vnl_vector type for this class. */
  typedef vnl_vector_fixed<TScalarType, itkGetStaticConstMacro(SpaceDimension)> InputVnlVectorType;
  typedef vnl_vector_fixed<TScalarType, itkGetStaticConstMacro(SpaceDimension)> OutputVnlVectorType;
  
  /** Standard coordinate point type for this class. */
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> InputPointType;
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> OutputPointType;

  /** Standard interpolated image. */
  typedef itk::Image< OutputPointType, NDimensions >  InterpolatedImageType;
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------



  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  /** This method sets the parameters for the transform value specified by the user. */
  void SetParameters(const ParametersType & parameters);

  /** This method set the spline paremeters so it is equivalent to an affine transform. */
  void SetAffineTransform(const ParametersType & parameters);

  /** Get the Transformation Parameters. */
  virtual const ParametersType& GetParameters(void) const;

  /** Get the pseudo-inverse transformation Parameters. */
  const ParametersType GetInverseParameters( void ) const;

  /** Solve a block-tridiagonal system by means of the Thomas algorithm */
  void Thomas( double* values, double* results, unsigned int* sizes, unsigned int order ) const;

  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------


  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  // Tranformation of a unique point:
  OutputPointType     TransformPoint( const InputPointType  &point ) const;
  // Transformation of vectors: purposely not implemented:
  virtual OutputVectorType TransformVector( const InputVectorType & ) const{
	  itkExceptionMacro( << "Method not aplicable for spline transform" );
  }
  virtual OutputVnlVectorType TransformVector( const InputVnlVectorType & ) const{
	  itkExceptionMacro( << "Method not aplicable for spline transform" );
  }
  virtual OutputCovariantVectorType TransformVector( const InputCovariantVectorType & ) const{
	  itkExceptionMacro( << "Method not aplicable for spline transform" );
  }
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  /** B-spline functions: */
  double Bi(double s, unsigned int i) const
  {
	  double res;
	  switch(i){
		  case 0:
			  res = (1.0-s)*(1.0-s)*(1.0-s)/6.0;
			  break;
		  case 1:
			  res = (3.0*s*s*s - 6.0*s*s + 4.0)/6.0;
			  break;
		  case 2:
			  res = (-3.0*s*s*s + 3.0*s*s + 3.0*s + 1.0)/6.0;
			  break;
		  case 3:
		  default:
			  res = s*s*s/6.0;
			  break;
	  }
	  return (res);
  }

  double B0(double s) const
  {
	  return ((1.0-s)*(1.0-s)*(1.0-s)/6.0);
  }

  double B1(double s) const 
  {
	  return ((3.0*s*s*s - 6.0*s*s + 4.0)/6.0);
  }

  double B2(double s) const
  {
	  return ((-3.0*s*s*s + 3.0*s*s + 3.0*s + 1.0)/6.0);
  }

  double B3(double s) const
  {
	  return (s*s*s/6.0);
  }
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------


  /** Compute the Jacobian Matrix of the transformation at one point */
  virtual const JacobianType & GetJacobian( const InputPointType  &point ) const;

  /** Compute the regularization penalty term */
  double GetCorrection();

  /** Compute the Jacobian of the regularization penalty term */
  typedef itk::Array<double> CorrectionJacobianType;
  CorrectionJacobianType GetCorrectionJacobian();

  /** Set the parameters to the IdentityTransform */
  void SetIdentity( void );

  /** Set the number of points in the grid */
  void SetNumberOfPoints( const unsigned int pointsx, const unsigned int pointsy, const unsigned int pointsz );

  /** Set the region the grid extends over */
  void SetGridRegion( InputPointType origin, InputVectorType spacing, SizeType size );

  /** Refine the grid to the next level of resolution */
  void RefineGrid( void );

  /** Interpolate given points */
  void InterpolateImage( InterpolatedImageType* image );

  /** Interpolate given points */
  typedef std::vector<InputPointType>  InputVector;
  typedef std::vector<OutputPointType> OutputVector;
  void InterpolatePoints( InputVector points, OutputVector values);
  
  /** Return the number of parameters that completely define the Transfom  */
  virtual unsigned int GetNumberOfParameters( void ) const;



protected:
  D3BSplineTransform();
  ~D3BSplineTransform();
  /** Print contents of an D3BSplineTransform. */
  void PrintSelf( std::ostream &os, Indent indent ) const;

private:
  D3BSplineTransform( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
  // PHI matrices that contains spline coefficients:
  double*** PHIx;
  double*** PHIy;
  double*** PHIz;
  // Parameters to return on the GetParameters() Method
  ParametersType* m_ReturnedParameters;
  // Number of points in each direction:
  unsigned int NX;
  unsigned int NY;
  unsigned int NZ;
  // Origin and extent of the grid:
  double m_Origin[3];
  double m_Extent[3];
  // Stores the upper-left corner of the last used 4 by 4 by 4 vicinity of the grid. It has to be mutable
  // because we need to change it from the const method GetJacobian()
  mutable unsigned int m_LastVisited[3];
  mutable JacobianType m_Jacobian;
}; //class D3BSplineTransform


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkD3BSplineTransform.txx"
#endif

#endif /* __itkD3BSplineTransform_h */
