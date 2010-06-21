/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkD1ForwardDifferencesFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkD1ForwardDifferencesFilter_h
#define __itkD1ForwardDifferencesFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT D1ForwardDifferencesFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef D1ForwardDifferencesFilter                      Self;
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
	
	typedef itk::Image< float, TInputImage::ImageDimension > DiffusionImageType;
	typedef typename DiffusionImageType::Pointer             DiffusionImagePointer;
	typedef typename DiffusionImageType::ConstPointer        DiffusionImageConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro( D1ForwardDifferencesFilter, ImageToImageFilter );

	/** Set and get the diffusion term field */
	void SetDiffusionTerm( DiffusionImagePointer term ){
		m_DiffusionTerm = term;
		return;
	}
	DiffusionImagePointer GetDiffusionTerm( void ){
		return( m_DiffusionTerm );
	}
	/** Set and get the time constant */
	itkSetMacro(                 TimeStep,   float        );
	itkGetConstReferenceMacro(   TimeStep,   float        );

	/** Set and get the filtering dimension */
	itkSetMacro(                 Direction,  unsigned int );
	itkGetConstReferenceMacro(   Direction,  unsigned int );
protected:
	D1ForwardDifferencesFilter();
	virtual ~D1ForwardDifferencesFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
private:
	D1ForwardDifferencesFilter(const Self&);          // Purposely not implemented
	void operator=(const Self&);              // Purposely not implemented
	float                 m_TimeStep;         // The time constant for anisotropic diffusion
	DiffusionImagePointer m_DiffusionTerm;    // The diffusion term field
	unsigned int          m_Direction;        // The dimension to filter
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkD1ForwardDifferencesFilter.txx"
#endif

#endif
