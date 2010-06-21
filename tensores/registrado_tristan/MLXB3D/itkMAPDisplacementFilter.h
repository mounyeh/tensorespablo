/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMAPDisplacementFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMAPDisplacementFilter_h
#define __itkMAPDisplacementFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template < class OutputPointType, unsigned int NDimension >
class ITK_EXPORT MAPDisplacementFilter : public ImageToImageFilter< itk::Image< itk::Array<double>, NDimension >, itk::Image< OutputPointType, NDimension > >
{
public:
	/** Convenient typedefs for simplifying declarations. */
	typedef itk::Image<itk::Array<double>, NDimension> InputImageType;
	typedef itk::Image<OutputPointType, NDimension>    OutputImageType;

	/** Standard class typedefs. */
	typedef MAPDisplacementFilter Self;
	typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
	typedef SmartPointer<Self> Pointer;
	typedef SmartPointer<const Self>  ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(MAPDisplacementFilter, ImageToImageFilter);

	/** Image typedef support. */
	typedef typename InputImageType::PixelType   InputPixelType;
	typedef typename OutputImageType::PixelType  OutputPixelType;
	
	typedef typename InputImageType::RegionType  InputImageRegionType;
	typedef typename OutputImageType::RegionType OutputImageRegionType;
	
	typedef typename InputImageType::SizeType    InputSizeType;
	typedef typename OutputImageType::SizeType   OutputSizeType;

	typedef typename InputImageType::SpacingType SpacingType;

	virtual void Initialize( SpacingType spacing, InputSizeType searchsize );

protected:
	MAPDisplacementFilter();
	void PrintSelf( std::ostream& os, Indent indent ) const;
	
	void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId );

private:
	MAPDisplacementFilter(const Self&); //purposely not implemented
	void operator=(const Self&);        //purposely not implemented
	itk::Array2D<double> m_Displacement;
	SpacingType          m_Spacing;
	InputSizeType        m_SearchSize;
	bool                 m_Init;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMAPDisplacementFilter.txx"
#endif

#endif
