/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeStrainScalars.h,v $
  Language:  C++
  Date:      $Date: 2008/02/7 14:28:51 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkComputeStrainScalars_h
#define __itkComputeStrainScalars_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
//#include "itkDWImages.h"

namespace itk
{
	
template < class TInputImage, class TOutputImage >
class ITK_EXPORT ComputeStrainScalars : public ImageToImageFilter< TInputImage, TOutputImage >
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
	typedef ComputeStrainScalars                                     Self;
	typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
	typedef SmartPointer<Self>                                    Pointer;
	typedef SmartPointer<const Self>                              ConstPointer;
	
	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	
	/** Run-time type information (and related methods). */
	itkTypeMacro( ComputeStrainScalars, ImageToImageFilter );
	
	/** Image typedef support. */
	typedef typename InputImageType::PixelType            InputPixelType;
	typedef typename OutputImageType::PixelType           OutputPixelType;
	typedef typename InputImageType::RegionType           InputImageRegionType;
	typedef typename OutputImageType::RegionType          OutputImageRegionType;
	typedef typename InputImageType::SizeType             InputSizeType;
	typedef typename InputPixelType::EigenValuesArrayType EigenValuesArrayType;
	
	/** The types of scalars to compute. */
	typedef enum ScalarParameter{INV,EIG0,EIG1,ST0,ST1,ST2} ScalarParameter;
//	typedef enum ScalarParameter{FA,RA,EIG0,EIG1,EIG2,MD,DC0,DC1,DC2,DC3,DC4,DC5,SC0,SC1,SC2} ScalarParameter;
	
	void SetComputeInvariant( void ){m_Scalar=INV;}

	void SetComputeEigVal( unsigned int ord )
	{
		switch( ord ){
			case 0:
				m_Scalar=EIG0;
				break;
			case 1:
				m_Scalar=EIG1;
				break;
			default:
				m_Scalar=EIG0;
				break;
		}
	}

	void SetComputeST( unsigned int ord )
	{
		switch( ord ){
			case 0:
				m_Scalar=ST0;
				break;
			case 1:
				m_Scalar=ST1;
				break;
			case 2:
				m_Scalar=ST2;
				break;
			default:
				m_Scalar=ST0;
				break;
		}
	}

	
protected:
	ComputeStrainScalars();
	virtual ~ComputeStrainScalars() {}
	// Threaded filter!
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
	ComputeStrainScalars(const Self&);   // purposely not implemented
	void operator=(const Self&);      // purposely not implemented
	// The scalar to compute:
	ScalarParameter m_Scalar;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeStrainScalars.txx"
#endif

#endif
