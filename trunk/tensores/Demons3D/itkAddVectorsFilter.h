/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAddVectorsFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAddVectorsFilter_h
#define __itkAddVectorsFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{
/** \class AddVectorsFilter
 * 
 */
template < class TScalar, unsigned int NDimension >
class ITK_EXPORT AddVectorsFilter :
	public ImageToImageFilter<  itk::Image< itk::Vector<TScalar,NDimension>, NDimension >, itk::Image< itk::Vector<TScalar,NDimension>, NDimension >  >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension,  unsigned int, NDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef itk::Image< itk::Vector<TScalar,NDimension>, NDimension >         InputImageType;
  typedef itk::Image< itk::Vector<TScalar,NDimension>, NDimension >         OutputImageType;

  /** Standard class typedefs. */
  typedef AddVectorsFilter                                                  Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >             Superclass;
  typedef SmartPointer<Self>                                                Pointer;
  typedef SmartPointer<const Self>                                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(  Self  );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AddVectorsFilter, ImageToImageFilter );
  
  /** Image typedef support. */
  typedef typename InputImageType::Pointer          InputImagePointer;
  typedef typename OutputImageType::Pointer         OutputImagePointer;

  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef typename OutputImageType::ConstPointer    OutputImageConstPointer;

  typedef typename InputImageType::PixelType        InputPixelType;
  typedef typename OutputImageType::PixelType       OutputPixelType;

  typedef typename InputImageType::RegionType       InputRegionType;
  typedef typename OutputImageType::RegionType      OutputRegionType;

  typedef typename InputImageType::SizeType         SizeType;
  typedef typename InputImageType::IndexType        IndexType;
  typedef typename InputImageType::SpacingType      SpacingType;

  /** Set the first input */
  void SetInput1( const InputImageType * image ){ this->SetInput( 0, image ); }
  
  /** Set the second input */
  void SetInput2( const InputImageType * image ){ this->SetInput( 1, image ); }

  /** Get the first input */
  const InputImageType * GetInput1( void ){ return this->GetInput( 0 ); }
  
  /** Get the second input */
  const InputImageType * GetInput2( void ){ return this->GetInput( 1 ); }

protected:
  AddVectorsFilter();
  virtual ~AddVectorsFilter() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;
  /** Multi-threaded filter */
  void ThreadedGenerateData( const OutputRegionType& outputRegionForThread, int threadId );
  /** We need a larger input than the output */
  virtual void GenerateInputRequestedRegion( ) throw( InvalidRequestedRegionError );

private:
  AddVectorsFilter(const Self&);   //purposely not implemented
  void operator=(const Self&);     //purposely not implemented
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAddVectorsFilter.txx"
#endif

#endif
