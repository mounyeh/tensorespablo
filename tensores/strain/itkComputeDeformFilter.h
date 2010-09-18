/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeDeformFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.16 $

=========================================================================*/
#ifndef __itkComputeDeformFilter_h
#define __itkComputeDeformFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"

namespace itk
{
	namespace Function {
		template< class TInput, class TOutput>
		class Deform
		{
		public:
			Deform() {}
			~Deform() {}
			inline TOutput operator()( const TInput & A )
			{
				return sqrt(A[0]*A[0]+A[1]*A[1]);
			}
		};
	}


template <class TInputImage, class TOutputImage>
class ITK_EXPORT ComputeDeformFilter :
	public UnaryFunctorImageFilter<  TInputImage,TOutputImage, Function::Deform< typename TInputImage::PixelType, typename TOutputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef ComputeDeformFilter           Self;
  typedef UnaryFunctorImageFilter<  TInputImage,  TOutputImage,  Function::Deform< typename TInputImage::PixelType, typename TOutputImage::PixelType>   > 
	                                Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
protected:
  ComputeDeformFilter() {}
  virtual ~ComputeDeformFilter() {}

private:
  ComputeDeformFilter(const Self&);      //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
