/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEigenValueCorrectionFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEigenValueCorrectionFilter_h
#define __itkEigenValueCorrectionFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT EigenValueCorrectionFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef EigenValueCorrectionFilter                       Self;
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

	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro(EigenValueCorrectionFilter, ImageToImageFilter);

	/** Set the inputs: */
	void SetInput1( const InputImageType* input ){
		this->SetInput( 0, input );
	}
	const InputImageType* GetInput1( void ){
		return this->GetInput( 0 );
	}
	void SetInput2( const InputImageType* input ){
		this->SetInput( 1, input );
	}
	const InputImageType* GetInput2( void ){
		return this->GetInput( 1 );
	}

	/** Set and get the beta parameter */
	itkSetMacro(                 Beta,   float   );
	itkGetConstReferenceMacro(   Beta,   float   );
protected:
	EigenValueCorrectionFilter();
	virtual ~EigenValueCorrectionFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
	EigenValueCorrectionFilter(const Self&); // purposely not implemented
	void operator=(const Self&);             // purposely not implemented
	float              m_Beta;               // The time constant for anisotropic diffusion
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEigenValueCorrectionFilter.txx"
#endif

#endif
