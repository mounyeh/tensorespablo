/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeLCCFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComputeLCCFilter_h
#define __itkComputeLCCFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{
/** \class ComputeLCCFilter
 * 
 */
template < const unsigned int NDimension >
class ITK_EXPORT ComputeLCCFilter :
	public ImageToImageFilter<    itk::Image< double, NDimension >, itk::Image< itk::Vector<double,NDimension>, NDimension >    >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension,  unsigned int, NDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef itk::Image< double, NDimension >                                InputImageType;
  typedef itk::Image< itk::Vector<double,NDimension>, NDimension >        OutputImageType;

  /** Standard class typedefs. */
  typedef ComputeLCCFilter                                                Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >           Superclass;
  typedef SmartPointer<Self>                                              Pointer;
  typedef SmartPointer<const Self>                                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(  Self  );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComputeLCCFilter, ImageToImageFilter );
  
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
  
  //-----------------------------------------------------------------------------------
  // Set and get the inputs:
  //    - Input 1: Fixed image ( I )
  //    - Input 2: Moving Image ( J )
  //    - Input 3: mean( I )
  //    - Input 4: mean( J )
  //    - Input 5: mean( I^2 )
  //    - Input 6: mean( J^2 )
  //    - Input 7: mean( I·J )
  /** Set the 1st input */
  void SetInput1(  const InputImageType * image ){ this->SetInput( 0, image ); }  
  /** Set the 2nd input */
  void SetInput2(  const InputImageType * image ){ this->SetInput( 1, image ); }
  /** Set the 3rd input */
  void SetInput3(  const InputImageType * image ){ this->SetInput( 2, image ); }
  /** Set the 4th input */
  void SetInput4(  const InputImageType * image ){ this->SetInput( 3, image ); }
  /** Set the 5th input */
  void SetInput5(  const InputImageType * image ){ this->SetInput( 4, image ); }
  /** Set the 6th input */
  void SetInput6(  const InputImageType * image ){ this->SetInput( 5, image ); }  
  /** Set the 7th input */
  void SetInput7(  const InputImageType * image ){ this->SetInput( 6, image ); }

  /** Get the 1st input */
  const InputImageType * GetInput1(  void ){ return this->GetInput( 0 ); }  
  /** Get the 2nd input */
  const InputImageType * GetInput2(  void ){ return this->GetInput( 1 ); }
  /** Get the 3rd input */
  const InputImageType * GetInput3(  void ){ return this->GetInput( 2 ); }
  /** Get the 4th input */
  const InputImageType * GetInput4(  void ){ return this->GetInput( 3 ); }
  /** Get the 5th input */
  const InputImageType * GetInput5(  void ){ return this->GetInput( 4 ); }
  /** Get the 6th input */
  const InputImageType * GetInput6(  void ){ return this->GetInput( 5 ); }  
  /** Get the 7th input */
  const InputImageType * GetInput7(  void ){ return this->GetInput( 6 ); }
  //-----------------------------------------------------------------------------------

  /** Get the global LCC value. Same as GetLCC */
  double GetMetricValue( void ) const{
	  return m_LCC;
  }

  /** Set and get tolerance for the local variances: */
  itkSetMacro(               Tol, double );
  itkGetConstReferenceMacro( Tol, double );

  /** Set and get tolerance for the global LCC: */
  itkSetMacro(               LCC, double  );
  itkGetConstReferenceMacro( LCC, double  );

protected:
  ComputeLCCFilter();
  virtual ~ComputeLCCFilter() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;
  /** Execute before filtering */
  void BeforeThreadedGenerateData();
  /** Multi-threaded filter */
  void ThreadedGenerateData( const OutputRegionType& outputRegionForThread, int threadId );
  /** Execute after filtering */
  void AfterThreadedGenerateData();
  /** We need a larger input than the output */
  virtual void GenerateInputRequestedRegion( ) throw( InvalidRequestedRegionError );

private:
  ComputeLCCFilter(const Self&);   //purposely not implemented
  void operator=(const Self&);     //purposely not implemented
  // The minimum variance in the vicinity:
  double               m_Tol;
  // Auxiliar value to sum the global LCC value:
  itk::Array< double > m_PerThreadLCC;
  // Global LCC value:
  double               m_LCC;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeLCCFilter.txx"
#endif

#endif
