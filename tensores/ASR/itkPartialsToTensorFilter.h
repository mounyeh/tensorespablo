/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPartialsToTensorFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPartialsToTensorFilter_h
#define __itkPartialsToTensorFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT PartialsToTensorFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef PartialsToTensorFilter                          Self;
	typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
	typedef SmartPointer< Self >                            Pointer;
	typedef SmartPointer< const Self >                      ConstPointer;
	
	/** Convenient typedefs for simplifying declarations. */
	typedef TInputImage                                     InputImageType;
	typedef typename InputImageType::Pointer                InputImagePointer;
	typedef typename InputImageType::ConstPointer           InputImageConstPointer;
	typedef TOutputImage                                    OutputImageType;
	typedef typename OutputImageType::Pointer               OutputImagePointer;
	
	/** Image typedef support. */
	typedef typename InputImageType::PixelType              InputPixelType;
	typedef typename OutputImageType::PixelType             OutputPixelType;
	typedef typename InputImageType::RegionType             InputImageRegionType;
	typedef typename OutputImageType::RegionType            OutputImageRegionType;
	typedef typename InputImageType::SizeType               InputSizeType;

	/** typedefs for the eigenvalue correction image */
	typedef itk::Image< float, TInputImage::ImageDimension >
		                                                    CorrectionType;
	typedef CorrectionType*                                 CorrectionPointer;
	typedef const CorrectionType*                           CorrectionConstPointer;
	
	void SetCorrection( CorrectionConstPointer correction ){
		m_Correction = const_cast< CorrectionPointer >(   correction   );
	}
	CorrectionPointer GetCorrection( void ){
		return m_Correction;
	}

	/** Set and get the correction image */


	/** Method for creation through the object factory. */
	itkNewMacro( Self );
	/** Run-time type information (and related methods). */
	itkTypeMacro( PartialsToTensorFilter, ImageToImageFilter );

	/** Set and get the time difference between eigenvalues */
	itkSetMacro(                 Difference,   float        );
	itkGetConstReferenceMacro(   Difference,   float        );

	/** Types for the diffusion tensor field */
	static const unsigned int TensorDimension = ( (TInputImage::ImageDimension)*(TInputImage::ImageDimension + 1) )/2;
protected:
	PartialsToTensorFilter();
	virtual ~PartialsToTensorFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
	PartialsToTensorFilter(const Self&); // purposely not implemented
	void operator=(const Self&);         // purposely not implemented
	float             m_Difference;      // Maximum difference between eigenvalues
	CorrectionPointer m_Correction;      // The correction for processed eigenvalues
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPartialsToTensorFilter.txx"
#endif

#endif
