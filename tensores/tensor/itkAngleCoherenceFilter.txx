/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAngleCoherenceFilter.txx,v $
  Language:  C++
  Date:      $Date: 2003/12/15 14:13:18 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAngleCoherenceFilter_txx
#define _itkAngleCoherenceFilter_txx

#include "itkAngleCoherenceFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"

namespace itk {


template<class TInputImage, class TMaskImage, class TOutputImage >
AngleCoherenceFilter<TInputImage, TMaskImage, TOutputImage>::AngleCoherenceFilter()
{
	// this filter requires two input images:
	this->SetNumberOfRequiredInputs(  1 );
	this->SetNumberOfRequiredOutputs( 1 );
	this->SetNumberOfThreads(1);
	
	m_BlockSize.Fill(  1 );
	m_Label=1;
}


template<class TInputImage, class TMaskImage, class TOutputImage >
void AngleCoherenceFilter<TInputImage, TMaskImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId ){
	/** ######################################################################################################### */
	// Counter:
	// Size of the array of doubles that compose the pixel:
	//unsigned long p_size = 1;
	unsigned long c_size = 1;	
	for( unsigned long p=0; p<TInputImage::ImageDimension; p++ ){
		//p_size *= (   2*( m_SearchSize[p] ) + 1   );
		c_size *= (   2*( m_BlockSize[p]  ) + 1   );
	}
	// Output pixel:
	OutputPixelType op;
	/** ######################################################################################################### */
	
	InputImageConstPointer input = this->GetInput();
	
	typedef itk::Image<itk::Vector<float,3>,3>	VectorImageType;
	VectorImageType::Pointer eigvectorImage = VectorImageType::New();

	eigvectorImage->CopyInformation( input );
	eigvectorImage->SetRegions( input->GetLargestPossibleRegion() );
	
	try{
		eigvectorImage->Allocate();
	}
	catch ( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
    }
	
	
	typedef itk::ImageRegionIterator<VectorImageType> VectorIteratorType;
	VectorIteratorType it_vec(eigvectorImage,outputRegionForThread);
	typedef itk::ImageRegionConstIterator< MaskImageType >       MaskIteratorType;
	MaskIteratorType   it_mask(this->GetMask(), outputRegionForThread);
	typedef itk::ImageRegionConstIteratorWithIndex<InputImageType>	InputIteratorType;
	InputIteratorType it1(input,outputRegionForThread);
	
	itk::SymmetricSecondRankTensor<float> Daux;
	
	itk::Vector<float,3> eigvec1;
	for (it1.GoToBegin(),it_vec.GoToBegin(),it_mask.GoToBegin();!it1.IsAtEnd();++it1, ++it_vec, ++it_mask) {
		if(it_mask.Get()==m_Label){
			itk::SymmetricSecondRankTensor<float>::EigenValuesArrayType eigval;
			itk::SymmetricSecondRankTensor<float>::EigenVectorsMatrixType eigvec;
			for(unsigned i=0; i<6; ++i){
				Daux[i]=it1.Get()[i];
			}
			Daux.ComputeEigenAnalysis( eigval, eigvec );
			eigvec1[0]=eigvec[2][0];
			eigvec1[1]=eigvec[2][1];
			eigvec1[2]=eigvec[2][2];
		} else {
			eigvec1.Fill(0);
		}
		it_vec.Set(eigvec1);
	}

	OutputImagePointer output      = this->GetOutput();
	
	// Iterators typedefs:
	typedef itk::ConstNeighborhoodIterator< VectorImageType >    VectorNeighIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutputIteratorType;

	// Boundary condition:
	ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;

	// Iterators creation:
	VectorNeighIteratorType itn( m_BlockSize, eigvectorImage, outputRegionForThread );
	OutputIteratorType it( output, outputRegionForThread );

	//itn.OverrideBoundaryCondition( &nbc );
	float coh;
	float cosangle;
	for ( it.GoToBegin(), itn.GoToBegin(), it_mask.GoToBegin(); !it.IsAtEnd(); ++it, ++it_mask, ++itn){
		if(it_mask.Get()!=m_Label){
			op = 0.0;
		}else{
			coh=0;
			for(unsigned long p=0; p<c_size; p++ ){
				// For all pixels in the block, calculate both mean values:
				cosangle = fabs(itn.GetCenterPixel()*itn.GetPixel(p));
				if(cosangle<1.0){
					coh += fabs(acos(cosangle));
				}else{
					coh += 0;
				}
			}
			op=coh/c_size;
		}
		it.Set( op );
		//--------------------------------------------------------------
	}
	//---------------------------------------------------------------------------------------------------------------
}


template<class TInputImage1, class TInputImage2, class TMaskImage>
void AngleCoherenceFilter<TInputImage1, TInputImage2, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf( os, indent );
	
	os << indent << "m_BlockSize: "   << m_BlockSize  << std::endl;
	
}


}// end namespace itk
#endif
