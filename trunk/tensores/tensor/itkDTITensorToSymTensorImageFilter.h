/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkDTItensorToSymTensorImageFilter_h
  Language:  C++
  Date:      2008/07/01
  Version:   
  
	 This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDTItensorToSymTensorImageFilter_h
#define __itkDTItensorToSymTensorImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{
  
/** \class itkDTItensorToSymmetricTensor
 * \brief Converts an DTITensor image into a symmetric second rank tensor image.
 * 
 *
 */
namespace Function {  
  
template< class TInput, class TOutput>
class DTITensorToSymTensor
{
public:
  typedef typename TInput::ComponentType      ComponentType;
  typedef typename itk::NumericTraits< ComponentType >::RealType  RealType;

  DTITensorToSymTensor() {}
  ~DTITensorToSymTensor() {}
  bool operator!=( const DTITensorToSymTensor & ) const
  {
    return false;
  }
  bool operator==( const DTITensorToSymTensor & other ) const
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
class ITK_EXPORT DTITensorToSymTensorImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Function::DTITensorToSymTensor< 
  typename TInputImage::PixelType, 
  typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef DTITensorToSymTensorImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Function::DTITensorToSymTensor< typename TInputImage::PixelType, 
                                                 typename TOutputImage::PixelType> >  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TInputImage::PixelType::ComponentType>));
  /** End concept checking */
#endif

protected:
  DTITensorToSymTensorImageFilter() {}
  virtual ~DTITensorToSymTensorImageFilter() {}

private:
  DTITensorToSymTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
