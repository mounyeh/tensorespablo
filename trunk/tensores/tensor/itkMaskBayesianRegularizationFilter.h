/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMaskBayesianRegularizationFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMaskBayesianRegularizationFilter_h
#define __itkMaskBayesianRegularizationFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{

template < class TMaskImage, const unsigned int NDimension >
class ITK_EXPORT MaskBayesianRegularizationFilter 
	: public ImageToImageFilter< itk::Image< itk::Array<double>, NDimension >,   itk::Image< itk::Array<double>, NDimension > >
{
public:
  /** Convenient typedefs for simplifying declarations. */
	typedef itk::Image<itk::Array<double>, NDimension> InputImageType;
	typedef itk::Image<itk::Array<double>, NDimension> OutputImageType;
	typedef TMaskImage								   MaskImageType;
	typedef typename TMaskImage::Pointer			   MaskImagePointer;

	/** Standard class typedefs. */
	typedef MaskBayesianRegularizationFilter                         Self;
	typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
	typedef SmartPointer<Self>                                   Pointer;
	typedef SmartPointer<const Self>                             ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(MaskBayesianRegularizationFilter, ImageToImageFilter);
	
	/** Image typedef support. */
	typedef typename InputImageType::PixelType   InputPixelType;
	typedef typename OutputImageType::PixelType  OutputPixelType;
	
	typedef typename InputImageType::RegionType  InputImageRegionType;
	typedef typename OutputImageType::RegionType OutputImageRegionType;
	
	typedef typename InputImageType::SizeType    InputSizeType;
	typedef typename OutputImageType::SizeType   OutputSizeType;

	/** MaskBayesianRegularizationFilter needs a larger input requested region than
	* the output requested region.  As such, MaskBayesianRegularizationFilter needs
	* to provide an implementation for GenerateInputRequestedRegion()
	* in order to inform the pipeline execution model.*/
	
	virtual void GenerateInputRequestedRegion( ) throw( InvalidRequestedRegionError );
	virtual void Initialize( );

	itkSetMacro(               SearchSize, InputSizeType  );
	itkGetConstReferenceMacro( SearchSize, InputSizeType  );

	itkSetMacro(               Sigma, itk::Array<double>  );
	itkGetConstReferenceMacro( Sigma, itk::Array<double>  );

	/** Set the mask input. */
	void SetMask( MaskImageType * image ){
		m_Mask = image;
	}

	/** Get mask input. */
	MaskImagePointer  GetMask( void ) const{
		return m_Mask;
	}
protected:
	MaskBayesianRegularizationFilter();
	virtual ~MaskBayesianRegularizationFilter() {}
	void PrintSelf(std::ostream& os, Indent indent) const;
	
	/** MaskBayesianRegularizationFilter can be implemented as a multithreaded filter.
	* Therefore, this implementation provides a ThreadedGenerateData()
	* routine which is called for each processing thread. The output
	* image data is allocated automatically by the superclass prior to
	* calling ThreadedGenerateData().  ThreadedGenerateData can only
	* write to the portion of the output image specified by the
	* parameter "outputRegionForThread"
	*
	* \sa ImageToImageFilter::ThreadedGenerateData(),
	*     ImageToImageFilter::GenerateData() */
	
	void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId );

private:
	MaskBayesianRegularizationFilter(const Self&); // Purposely not implemented
	void operator=(const Self&);               // Purposely not implemented
	InputSizeType        m_Radius;             // Size of the regularization vicinity: 3x3x...x3, so m_Radius=[1, 1, ..., 1]
	InputSizeType        m_SearchSize;         // Radius of the block matching search vicinity
	itk::Array2D<double> m_Distance;           // Euclidean distances among each pair of displacement vectors
	itk::Array<double>   m_Sigma;              // Variance of the algorithm 
	bool                 m_Init;               // Indicator of algorithm initiallization
	MaskImagePointer	 m_Mask;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMaskBayesianRegularizationFilter.txx"
#endif

#endif
