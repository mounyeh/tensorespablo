/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkModifiedLinearInterpolateImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2003/11/16 04:38:10 $
  Version:   $Revision: 1.27 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkModifiedLinearInterpolateImageFunction_h
#define _itkModifiedLinearInterpolateImageFunction_h

#include "itkLinearInterpolateImageFunction.h"
#include "itkPoint.h"

namespace itk
{

/** \class ModifiedLinearInterpolateImageFunction
 * \brief Linearly interpolate an image at specified positions.
 *
 * ModifiedLinearInterpolateImageFunction linearly interpolates image intensity at
 * a non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type 
 * (e.g. float or double).
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for images with scalar pixel
 * types. For vector images use VectorModifiedLinearInterpolateImageFunction.
 *
 * \sa VectorModifiedLinearInterpolateImageFunction
 *
 * \ingroup ImageFunctions ImageInterpolators 
 */
template <class TInputImage, class TCoordRep = float>
class ITK_EXPORT ModifiedLinearInterpolateImageFunction : 
  public LinearInterpolateImageFunction<TInputImage,TCoordRep> 
{
public:
  /** Standard class typedefs. */
  typedef ModifiedLinearInterpolateImageFunction                  Self;
  typedef LinearInterpolateImageFunction<TInputImage,TCoordRep>   Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(ModifiedLinearInterpolateImageFunction, InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType                         OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType                     InputImageType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType                           RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType                          IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType                ContinuousIndexType;

  typedef Point<TCoordRep,itkGetStaticConstMacro(ImageDimension)> PointType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a 
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( const ContinuousIndexType & index ) const;

  /** Overload IsInsideBuffer() to achieve the desired behaviour when the point is outside the image: */
  bool IsInsideBuffer( const PointType & point ) const;

protected:
  ModifiedLinearInterpolateImageFunction();
  ~ModifiedLinearInterpolateImageFunction(){};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  ModifiedLinearInterpolateImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;  

};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkModifiedLinearInterpolateImageFunction.txx"
#endif

#endif
