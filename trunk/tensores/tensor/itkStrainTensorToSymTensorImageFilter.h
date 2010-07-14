/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkStraintensorToSymTensorImageFilter_h
  Language:  C++
  Date:      2008/07/01
  Version:   
  
	 This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkStraintensorToSymTensorImageFilter_h
#define __itkStraintensorToSymTensorImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{
  
/** \class itkStraintensorToSymmetricTensor
 * \brief Converts a StrainTensor image into a symmetric second rank tensor image.
 * 
 *
 */
namespace Function {  
  
template< class TInput, class TOutput>
class StrainTensorToSymTensor
{
public:
//  typedef typename TInput::ComponentType      ComponentType;
  typedef typename itk::NumericTraits< float >::RealType  RealType;

  StrainTensorToSymTensor() {}
  ~StrainTensorToSymTensor() {}
  bool operator!=( const StrainTensorToSymTensor & ) const
  {
    return false;
  }
  bool operator==( const StrainTensorToSymTensor & other ) const
  {
    return !(*this != other);
  }
  inline TOutput operator()( const TInput & A )
  {	
	TOutput outputData;
	for(unsigned i=0; i<A.GetNumberOfComponents(); i++){
		outputData[i]=A[i];
		}
   return outputData; }
}; 
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT StrainTensorToSymTensorImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Function::StrainTensorToSymTensor< 
  typename TInputImage::PixelType, 
  typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef StrainTensorToSymTensorImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Function::StrainTensorToSymTensor< typename TInputImage::PixelType, 
                                                 typename TOutputImage::PixelType> >  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<float>));
  /** End concept checking */
#endif

protected:
  StrainTensorToSymTensorImageFilter() {}
  virtual ~StrainTensorToSymTensorImageFilter() {}

private:
  StrainTensorToSymTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
