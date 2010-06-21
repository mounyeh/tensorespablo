/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkModifiedLinearInterpolateImageFunction.txx,v $
  Language:  C++
  Date:      $Date: 2004/12/12 22:07:23 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkModifiedLinearInterpolateImageFunction_txx
#define _itkModifiedLinearInterpolateImageFunction_txx
#include "itkModifiedLinearInterpolateImageFunction.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Define the number of neighbors
 */
template<class TInputImage, class TCoordRep>
const unsigned long
ModifiedLinearInterpolateImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep>
ModifiedLinearInterpolateImageFunction< TInputImage, TCoordRep >
::ModifiedLinearInterpolateImageFunction()
{

}

template<class TInputImage, class TCoordRep>
bool ModifiedLinearInterpolateImageFunction< TInputImage, TCoordRep >
::IsInsideBuffer( const PointType & point ) const
{
	return true;
}


/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep>
void
ModifiedLinearInterpolateImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename ModifiedLinearInterpolateImageFunction< TInputImage, TCoordRep >::OutputType
ModifiedLinearInterpolateImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex( const ContinuousIndexType& index ) const
{
	// The point is inside the buffer: default behaviour:
	if( Superclass::IsInsideBuffer(index) ){
		return Superclass::EvaluateAtContinuousIndex( index );
	}
	// The point is outside the image:
	ContinuousIndexType extremeIndex = index;
	for( unsigned int dim = 0; dim < ImageDimension; dim++ ){
		if( index[dim]<Superclass::GetStartContinuousIndex()[dim] )
			extremeIndex[dim] = Superclass::GetStartContinuousIndex()[dim];
		if( index[dim]>Superclass::GetEndContinuousIndex()[dim] )
			extremeIndex[dim] = Superclass::GetEndContinuousIndex()[dim];
	}
	return Superclass::EvaluateAtContinuousIndex( extremeIndex );
}

} // namespace itk

#endif
