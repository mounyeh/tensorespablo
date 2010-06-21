/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSRADFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSRADFilter_h
#define __itkSRADFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkSRADNDThomasFilter.h"
#include "itkSRADDiffusionTensorFilter.h"
#include "itkExpImageFilter.h"
#include "itkLogImageFilter.h"
#include "itkScaleGrayLevelFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT SRADFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef SRADFilter                                      Self;
	typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
	typedef SmartPointer< Self >                            Pointer;
	typedef SmartPointer< const Self >                      ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );
	/** Run-time type information (and related methods). */
	itkTypeMacro( SRADFilter, ImageToImageFilter );
	
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

	/** Types for the diffusion tensor field */
	typedef itk::Image< float,TInputImage::ImageDimension > TensorImageType;
	typedef typename TensorImageType::Pointer               TensorImagePointer;
	typedef typename TensorImageType::ConstPointer          TensorImageConstPointer;
	
	/** Types for the auxiliar filters */
	typedef itk::ExpImageFilter< InternalImageType, InternalImageType >            ExpType;
	typedef typename ExpType::Pointer                                              ExpPointer;
	typedef itk::LogImageFilter< InternalImageType, InternalImageType >            LogType;
	typedef typename LogType::Pointer                                              LogPointer;
	typedef itk::ScaleGrayLevelFilter< InputImageType, InternalImageType>          InScaleType;
	typedef typename InScaleType::Pointer                                          InScalePointer;
	typedef itk::ScaleGrayLevelFilter< InternalImageType, OutputImageType >        OutScaleType;
	typedef typename OutScaleType::Pointer                                         OutScalePointer;
	typedef itk::SRADNDThomasFilter< InternalImageType, InternalImageType >        ThomasType;
	typedef typename ThomasType::Pointer                                           ThomasPointer;
	typedef itk::SRADDiffusionTensorFilter< InternalImageType,TensorImageType > TensorFilterType;
	typedef typename TensorFilterType::Pointer                                     TensorFilterPointer;

	/** Set and get the number of iterations */
	itkSetMacro(               Iter, unsigned int );
	itkGetConstReferenceMacro( Iter, unsigned int );

	/** Set and get the time constant */
	itkSetMacro(               TimeStep, float );
	itkGetConstReferenceMacro( TimeStep, float );

	/** Set and get the time constant */
	itkSetMacro(               ExpConstant, float );
	itkGetConstReferenceMacro( ExpConstant, float );
	
	/** Set and get the stopping condition parameter */
	itkSetMacro(                 Beta,       float          );
	itkGetConstReferenceMacro(   Beta,       float          );

	itkSetMacro(                 Q0,         double         );
	itkGetConstReferenceMacro(   Q0,         double         );
	
	/** Set and get the diffusion tensor filter */
	void SetTensorFilter( TensorFilterType* tensor ){
		m_Tensor = const_cast< TensorFilterType* >( tensor );
	}
	TensorFilterType* GetTensorFilter( void ){
		return m_Tensor;
	}
protected:
	SRADFilter();
	~SRADFilter();
	void PrintSelf( std::ostream& os, Indent indent ) const;
	void GenerateInputRequestedRegion();
	void EnlargeOutputRequestedRegion( DataObject *output );
	void GenerateData( );
private:
	SRADFilter(const Self&);           // Purposely not implemented
	void operator=(const Self&);       // Purposely not implemented
	unsigned int        m_Iter;        // The number of iterations to perform
	float               m_TimeStep;    // Time step for the diffusion process
	float               m_ExpConstant; // Log-compound constant
	TensorFilterPointer m_Tensor;      // The filter to compute the diffusion tensor
	float  m_Beta;                     // Eigenvalue correction parameter
	double m_Q0;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSRADFilter.txx"
#endif

#endif

