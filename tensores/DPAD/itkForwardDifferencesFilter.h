/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkForwardDifferencesFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkForwardDifferencesFilter_h
#define __itkForwardDifferencesFilter_h

#include "itkImageToImageFilter.h"
#include "itkD1ForwardDifferencesFilter.h"
#include "itkAddImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ForwardDifferencesFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef ForwardDifferencesFilter                        Self;
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
	itkTypeMacro( ForwardDifferencesFilter, ImageToImageFilter );

	/** Types for the diffusion term */
	typedef itk::Image< float, TInputImage::ImageDimension > DiffusionImageType;
	typedef typename DiffusionImageType::Pointer             DiffusionImagePointer;
	typedef typename DiffusionImageType::ConstPointer        DiffusionImageConstPointer;

	/** Types for the auxiliar filters */
	typedef itk::D1ForwardDifferencesFilter< TInputImage, TOutputImage >
		                                                    ForwardType;
	typedef typename ForwardType::Pointer                   ForwardPointer;
	typedef itk::AddImageFilter< TOutputImage, TOutputImage, TOutputImage >
		                                                    AddType;
	typedef typename AddType::Pointer                       AddPointer;

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
protected:
	ForwardDifferencesFilter();
	virtual ~ForwardDifferencesFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void GenerateData( );
	void EnlargeOutputRequestedRegion(DataObject *output);
private:
	ForwardDifferencesFilter(const Self&);  // Purposely not implemented
	void operator=(const Self&);            // Purposely not implemented
	float                 m_TimeStep;       // The time constant for anisotropic diffusion
	DiffusionImagePointer m_DiffusionTerm;  // The diffusion term field
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkForwardDifferencesFilter.txx"
#endif

#endif
