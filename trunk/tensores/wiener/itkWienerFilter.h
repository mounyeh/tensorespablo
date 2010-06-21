/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDWTFilter.h $
  Language:  C++
  Date:      $Date: 2006/04/25 $
  Version:   $Revision: 1.0 $

=========================================================================*/
#ifndef __itkWienerFilter_h
#define __itkWienerFilter_h

#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{

template <class TInputImage,class TOutputImage>
class ITK_EXPORT WienerFilter :
	public ImageToImageFilter< TInputImage, TOutputImage > 
{
public:
  /** Standard class typedefs. */
  typedef WienerFilter                                      Self;
  typedef ImageToImageFilter< TInputImage,TOutputImage >    Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Type macro that defines a name for this class */
  itkTypeMacro( WienerFilter, ImageToImageFilter );

  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Type of the input and output image */
  typedef TInputImage                                       InputImageType;
  typedef typename InputImageType::Pointer                  InputImagePointer;
  typedef typename InputImageType::ConstPointer             InputImageConstPointer;
  typedef typename InputImageType::RegionType               InputImageRegionType;
  typedef typename InputImageType::SizeType                 InputImageSizeType;
  typedef typename InputImageType::IndexType                InputImageIndexType;

  typedef TOutputImage                                      OutputImageType;
  typedef typename OutputImageType::Pointer                 OutputImagePointer;
  typedef typename OutputImageType::RegionType              OutputImageRegionType;
  typedef typename OutputImageType::SizeType                OutputImageSizeType;
  typedef typename OutputImageType::IndexType               OutputImageIndexType;
  typedef typename OutputImageType::PixelType               OutputPixelType;
  
  typedef typename itk::Image<float,InputImageDimension> FloatImageType; 
  typedef typename FloatImageType::Pointer FloatImagePointer; 

  typedef typename itk::ConstNeighborhoodIterator<FloatImageType> NeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::RadiusType myRadiusType; 
  /** Unidimensional wavelet filter type */
  void SetInput( const TInputImage * image);
  itkSetMacro(opcion, unsigned int);
  itkGetConstReferenceMacro(opcion, unsigned int );
 /** Set and get lambda : */
  itkSetMacro(               lambda,  float );
  itkGetConstReferenceMacro( lambda,  float );
  /** Set and get the filter iterations: */
  itkSetMacro(               iterations, unsigned int );
  itkGetConstReferenceMacro( iterations, unsigned int );
    
  itkSetMacro(               gamma,  float );
  itkGetConstReferenceMacro( gamma,  float );

 
  itkSetMacro( Noise, bool );  
  itkSetMacro( Bias, bool );


protected:
  WienerFilter();
  virtual ~WienerFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  //void GenerateOutputInformation( void );
  //void GenerateInputRequestedRegion( ) throw( InvalidRequestedRegionError );
  //virtual void GenerateOutputRequestedRegion(DataObject *output);
  void GenerateData( );

private:  
  WienerFilter(const Self&);      //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  float              m_lambda;
  unsigned int       m_iterations;
  float m_gamma;
  unsigned int m_opcion;
  bool m_Noise;
  bool m_Bias;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWienerFilter.txx"
#endif

#endif

