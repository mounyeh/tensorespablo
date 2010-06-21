/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDPADFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDPADFilter_h
#define __itkDPADFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkDPADNDThomasFilter.h"
#include "itkForwardDifferencesFilter.h"
#include "itkComputeDiffusionTermFilter.h"
#include "itkCastImageFilter.h"


namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT DPADFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef DPADFilter                                      Self;
	typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
	typedef SmartPointer< Self >                            Pointer;
	typedef SmartPointer< const Self >                      ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );
	/** Run-time type information (and related methods). */
	itkTypeMacro( DPADFilter, ImageToImageFilter );
	
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

	/** Internal Image Type */
	typedef itk::Image< float,TInputImage::ImageDimension > InternalImageType;
	
	/** Types for the auxiliar filters */
	typedef itk::DPADNDThomasFilter< InternalImageType, InternalImageType >        ThomasType;
	typedef typename ThomasType::Pointer                                           ThomasPointer;
	typedef itk::ForwardDifferencesFilter< InternalImageType, InternalImageType >  ForwardDifferencesType;
	typedef typename ForwardDifferencesType::Pointer                               ForwardDifferencesPointer;
	typedef itk::ComputeDiffusionTermFilter< InputImageType, InternalImageType >   DiffusionTermType;
	typedef typename DiffusionTermType::Pointer                                    DiffusionTermPointer;
	typedef itk::CastImageFilter<InputImageType, InternalImageType>                InCastType;
	typedef typename InCastType::Pointer                                           InCastPointer;
	typedef itk::CastImageFilter<InternalImageType, OutputImageType>               OutCastType;
	typedef typename OutCastType::Pointer                                          OutCastPointer;

	/** Set and get the number of iterations */
	itkSetMacro(               Iter, unsigned int );
	itkGetConstReferenceMacro( Iter, unsigned int );

	/** Set and get the time constant */
	itkSetMacro(               TimeStep, float );
	itkGetConstReferenceMacro( TimeStep, float );
	
	/** Set and get the radius of the vicinity */
	itkSetMacro(               Radius, InputSizeType );
	itkGetConstReferenceMacro( Radius, InputSizeType );
	
	/** Set and get the diffusion term filter */
	void SetDiffusionTerm( DiffusionTermType* term ){
		m_DiffusionTerm = const_cast< DiffusionTermType* >( term );
	}
	DiffusionTermType* GetDiffusionTerm( void ){
		return m_DiffusionTerm;
	}
	itkSetMacro(                 UseMode,    bool  );
	itkGetConstReferenceMacro(   UseMode,    bool  );
	void SetUseModeOn(){
		this->SetUseMode( true );
	}
	void SetUseModeOff(){
		this->SetUseMode( false );
	}
	itkSetMacro(                 UseAOS,    bool  );
	itkGetConstReferenceMacro(   UseAOS,    bool  );
	void SetUseAOSOn(){
		this->SetUseAOS( true );
	}
	void SetUseAOSOff(){
		this->SetUseAOS( false );
	}
	void SetUseModeDiffusionTerm( DiffusionTermMode mode )
	{
		m_DiffusionTermMode = mode;
	}
	DiffusionTermMode GetUseModeDiffusionTerm( void ){
		return m_DiffusionTermMode;
	}	
	void SetUseAjaDiffusionTerm( void )
	{
		m_DiffusionTermMode = __Aja__;
	}
	void SetUseYuDiffusionTerm( void )
	{
		m_DiffusionTermMode = __Yu__;
	}
	void SetUseSimplifiedDiffusionTerm( void )
	{
		m_DiffusionTermMode = __Simplified__;
	}
protected:
	DPADFilter();
	~DPADFilter();
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void EnlargeOutputRequestedRegion( DataObject *output );
	void GenerateData( );
private:
	DPADFilter(const Self&);                 // Purposely not implemented
	void operator=(const Self&);             // Purposely not implemented
	unsigned int          m_Iter;            // The number of iterations to perform
	float                 m_TimeStep;        // Time step for the diffusion process
	DiffusionTermPointer  m_DiffusionTerm;   // The filter to compute the diffusion term
	InputSizeType         m_Radius;          // Vicinity radius (for future use in subclasses)
	bool                  m_UseMode;
	bool                  m_UseAOS;
	DiffusionTermMode     m_DiffusionTermMode;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDPADFilter.txx"
#endif

#endif

