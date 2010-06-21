/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDPADNDThomasFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDPADNDThomasFilter_h
#define __itkDPADNDThomasFilter_h

#include "itkImageToImageFilter.h"
#include "itkDPAD1DThomasFilter.h"
#include "itkAddImageFilter.h"
#include "itkImage.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT DPADNDThomasFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef DPADNDThomasFilter                              Self;
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
	itkTypeMacro( DPADNDThomasFilter, ImageToImageFilter );

	typedef itk::Image< float, TInputImage::ImageDimension > DiffusionImageType;
	typedef typename DiffusionImageType::Pointer             DiffusionImagePointer;
	typedef typename DiffusionImageType::ConstPointer        DiffusionImageConstPointer;
	
	/** Types for the auxiliar filters */
	typedef itk::DPAD1DThomasFilter< TInputImage, TOutputImage >
		                                                    ThomasType;
	typedef typename ThomasType::Pointer                    ThomasPointer;
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
	DPADNDThomasFilter();
	virtual ~DPADNDThomasFilter() {}
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void EnlargeOutputRequestedRegion(DataObject *output);
	void GenerateData( );
private:
	DPADNDThomasFilter(const Self&);        // Purposely not implemented
	void operator=(const Self&);            // Purposely not implemented
	float                  m_TimeStep;      // The time constant for anisotropic diffusion
	DiffusionImagePointer  m_DiffusionTerm; // The diffusion term field
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDPADNDThomasFilter.txx"
#endif

#endif
