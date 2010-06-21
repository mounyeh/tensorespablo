/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeDTIScalars.h,v $
  Language:  C++
  Date:      $Date: 2008/02/7 14:28:51 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkComputeDTIScalars_h
#define __itkComputeDTIScalars_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkDWImages.h"

namespace itk
{
	
template < class TInputImage, class TOutputImage >
class ITK_EXPORT ComputeDTIScalars : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Convenient typedefs for simplifying declarations. */
	typedef TInputImage                                       InputImageType;
	typedef typename InputImageType::Pointer                  InputImagePointer;
	typedef typename InputImageType::ConstPointer             InputImageConstPointer;
	typedef TOutputImage                                      OutputImageType;
	typedef typename OutputImageType::Pointer                 OutputImagePointer;
	typedef typename OutputImageType::ConstPointer            OutputImageConstPointer;
	
	/** Standard class typedefs. */
	typedef ComputeDTIScalars                                     Self;
	typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
	typedef SmartPointer<Self>                                    Pointer;
	typedef SmartPointer<const Self>                              ConstPointer;
	
	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	
	/** Run-time type information (and related methods). */
	itkTypeMacro( ComputeDTIScalars, ImageToImageFilter );
	
	/** Image typedef support. */
	typedef typename InputImageType::PixelType            InputPixelType;
	typedef typename OutputImageType::PixelType           OutputPixelType;
	typedef typename InputImageType::RegionType           InputImageRegionType;
	typedef typename OutputImageType::RegionType          OutputImageRegionType;
	typedef typename InputImageType::SizeType             InputSizeType;
	typedef typename InputPixelType::EigenValuesArrayType EigenValuesArrayType;
	
	/** The types of scalars to compute. */
	typedef enum ScalarParameter{FA,RA,EIG0,EIG1,EIG2,MD,DC0,DC1,DC2,DC3,DC4,DC5,SC0,SC1,SC2} ScalarParameter;
	
	void SetComputeFA( void ){m_Scalar=FA;}
	void SetComputeRA( void ){m_Scalar=RA;}
	void SetComputeEigVal( unsigned int ord )
	{
		switch( ord ){
			case 0:
				m_Scalar=EIG0;
				break;
			case 1:
				m_Scalar=EIG1;
				break;
			case 2:
			default:
				m_Scalar=EIG2;
				break;
		}
	}
	void SetComputeMD( void ){m_Scalar=MD;}
	void SetComputeDC( unsigned int ord )
	{
		switch( ord ){
			case 0:
				m_Scalar=DC0;
				break;
			case 1:
				m_Scalar=DC1;
				break;
			case 2:
				m_Scalar=DC2;
				break;
			case 3:
				m_Scalar=DC3;
				break;
			case 4:
				m_Scalar=DC4;
				break;
			case 5:
			default:
				m_Scalar=DC5;
				break;
		}
	}
	void SetComputeShapeCoefficients( unsigned int ord )
	{
		switch( ord ){
			case 0:
				m_Scalar=SC0;
				break;
			case 1:
				m_Scalar=SC1;
				break;
			case 2:
			default:
				m_Scalar=SC2;
				break;
		}
	}
protected:
	ComputeDTIScalars();
	virtual ~ComputeDTIScalars() {}
	// Threaded filter!
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
	ComputeDTIScalars(const Self&);   // purposely not implemented
	void operator=(const Self&);      // purposely not implemented
	// The scalar to compute:
	ScalarParameter m_Scalar;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeDTIScalars.txx"
#endif

#endif
