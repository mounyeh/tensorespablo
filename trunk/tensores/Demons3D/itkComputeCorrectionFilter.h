/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeCorrectionFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComputeCorrectionFilter_h
#define __itkComputeCorrectionFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{
/** \class ComputeCorrectionFilter
 * 
 */
template < unsigned int NDimension >
class ITK_EXPORT ComputeCorrectionFilter :
	public ImageToImageFilter<  itk::Image< itk::Vector<double,NDimension>, NDimension >, itk::Image< itk::Vector<double,NDimension>, NDimension >  >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension,  unsigned int, NDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef itk::Image< itk::Vector<double,NDimension>, NDimension >          InputImage1Type;
  typedef itk::Image< itk::Vector<double,NDimension>, NDimension >          InputImage2Type;
  typedef itk::Image< itk::Vector<double,NDimension>, NDimension >          OutputImageType;

  /** Standard class typedefs. */
  typedef ComputeCorrectionFilter                                           Self;
  typedef ImageToImageFilter< InputImage1Type, OutputImageType >            Superclass;
  typedef SmartPointer<Self>                                                Pointer;
  typedef SmartPointer<const Self>                                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(  Self  );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComputeCorrectionFilter, ImageToImageFilter );
  
  /** Image typedef support. */
  typedef typename InputImage1Type::Pointer         InputImage1Pointer;
  typedef typename InputImage1Type::Pointer         InputImage2Pointer;
  typedef typename OutputImageType::Pointer         OutputImagePointer;

  typedef typename InputImage1Type::ConstPointer    InputImage1ConstPointer;
  typedef typename InputImage1Type::ConstPointer    InputImage2ConstPointer;
  typedef typename OutputImageType::ConstPointer    OutputImageConstPointer;

  typedef typename InputImage1Type::PixelType       Input1PixelType;
  typedef typename InputImage2Type::PixelType       Input2PixelType;
  typedef typename OutputImageType::PixelType       OutputPixelType;

  typedef typename InputImage1Type::RegionType      Input1RegionType;
  typedef typename InputImage2Type::RegionType      Input2RegionType;
  typedef typename OutputImageType::RegionType      OutputRegionType;

  typedef typename InputImage1Type::SizeType        SizeType;
  typedef typename InputImage2Type::SizeType        SizeType2;

  typedef typename InputImage1Type::IndexType       IndexType;
  typedef typename InputImage2Type::IndexType       IndexType2;

  typedef typename InputImage1Type::SpacingType     SpacingType;

  /** Set the first input */
  void SetInput1( const InputImage1Type * image ){ this->SetInput( 0, image ); }
  
  /** Set the second input */
  void SetInput2( const InputImage2Type * image ){ this->SetInput( 1, image ); }

  /** Get the first input */
  const InputImage1Type * GetInput1( void ){ return this->GetInput( 0 ); }
  
  /** Get the second input */
  const InputImage2Type * GetInput2( void ){ return this->GetInput( 1 ); }

  /** Set and get the Levenberg-Marquardt parameter: */
  itkSetMacro(               Lambda, double );
  itkGetConstReferenceMacro( Lambda, double );

protected:
  ComputeCorrectionFilter();
  virtual ~ComputeCorrectionFilter() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;
  /** Execute before filtering */
  void BeforeThreadedGenerateData();
  /** Multi-threaded filter */
  void ThreadedGenerateData( const OutputRegionType& outputRegionForThread, int threadId );
  /** We need a larger input than the output */
  virtual void GenerateInputRequestedRegion( ) throw( InvalidRequestedRegionError );

private:
  ComputeCorrectionFilter(const Self&);   //purposely not implemented
  void operator=(const Self&);            //purposely not implemented
  double      m_Lambda;
  SpacingType m_Spacing;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeCorrectionFilter.txx"
#endif

#endif
