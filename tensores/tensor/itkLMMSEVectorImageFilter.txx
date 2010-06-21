/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLMMSEVectorImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.14 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkLMMSEVectorImageFilter_txx
#define _itkLMMSEVectorImageFilter_txx
#include "itkLMMSEVectorImageFilter.h"
#ifndef M_PI
#define M_PI 3.141516
#endif

namespace itk
{

template <class TInputImage, class TOutputImage>
LMMSEVectorImageFilter<TInputImage, TOutputImage>
::LMMSEVectorImageFilter()
{

  m_RadiusEstimation.Fill( 1 );
  m_RadiusFiltering.Fill( 1 );
  m_Iterations = 5;
  //m_MinimumNumberOfUsedVoxelsEstimation = 1;
  m_MinimumNumberOfUsedVoxelsFiltering = 1;
  m_FirstBaseline = 0;
  m_UseAbsoluteValue = false;
  m_KeepValue = false;
  m_MinimumNoiseSTD = 1;
  m_MaximumNoiseSTD = 10000;
  //m_HistogramResolutionFactor = 2.0;
  //m_MaximumNumberOfBins = 200000;
}


template <class TInputImage, class TOutputImage>
void LMMSEVectorImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  // Get pointer to the input
  InputImagePointer inputPtr = const_cast< TInputImage * >( this->GetInput() );
  inputPtr->SetRequestedRegion( inputPtr->GetLargestPossibleRegion() );
}

template <class TInputImage, class TOutputImage>
void LMMSEVectorImageFilter< TInputImage, TOutputImage>
::GenerateData()
{
	// Auxiliar filters:
	LocalStatsPointer stats      = LocalStatsType::New();
	ThresholdPointer  threshold  = ThresholdType::New();
	HistogramPointer  histogram  = HistogramType::New();
	LMMSEStepPointer  lmmse      = LMMSEStepType::New();
	// For each iteration:
	for( unsigned int iter=0; iter<m_Iterations; ++iter ){
		// Create a minipieline to compute the noise:
		if( iter==0 )
			stats->SetInput( this->GetInput() );
		else
			stats->SetInput( lmmse->GetOutput() );
		typename LocalStatsType::IndicatorType indicator( 1 );
		// Compute the local characteristic:
		indicator[0] = m_FirstBaseline;
		stats->SetRadius( m_RadiusEstimation );
		stats->SetChannels( m_Channels );
		stats->SetUseNeighborhoodBaselines();
		stats->SetIndicator( indicator );
		stats->Update();
		// Compute the threshold in order to maximize the inter-class separability:
		threshold->SetMin( stats->GetMin() );
		threshold->SetMax( stats->GetMax() );
		threshold->SetW( 2.0f );
		threshold->SetBins( 2048 );
		threshold->SetInput( stats->GetOutput() );
		threshold->Update();
		double th = threshold->GetThreshold();
		// Create the minipipeline to compute the histogram:
		histogram->SetInput( stats->GetOutput() );
		histogram->SetMin(  2.0f );
		histogram->SetMax(   th  );
		histogram->SetBins( 256 );
		histogram->Update();
		// Compute the noise:
		double noiseStd  = histogram->GetMode();
		noiseStd        *= sqrt(2/M_PI);
std::cerr << "###########################################################" << std::endl;
std::cerr << "Estimated noise is: sigma = " << noiseStd << std::endl;
std::cerr << "###########################################################" << std::endl;
		if( noiseStd>m_MaximumNoiseSTD )
			noiseStd = m_MaximumNoiseSTD;
		if( noiseStd<m_MinimumNoiseSTD )
			noiseStd = m_MinimumNoiseSTD;
		// Now, set the LMMSE estimator input:
		if( iter==0 )
			lmmse->SetInput(  this->GetInput()  );
		else
			lmmse->SetInput(  lmmse->GetOutput()  );
		// Set the parameters:
		lmmse->SetFirstBaseline( m_FirstBaseline );
		lmmse->SetUseAbsoluteValue( m_UseAbsoluteValue );
		lmmse->SetKeepValue( m_KeepValue );
		lmmse->SetMinimumNumberOfUsedVoxelsFiltering( m_MinimumNumberOfUsedVoxelsFiltering );
		lmmse->SetNoiseVariance( noiseStd*noiseStd );
		lmmse->SetRadius( m_RadiusFiltering );
		lmmse->SetChannels( m_Channels );
		lmmse->Modified();
		lmmse->Update();
	}
	this->GraftOutput( lmmse->GetOutput() );
	return;
}


/** Standard "PrintSelf" method */
template <class TInputImage, class TOutput>
void LMMSEVectorImageFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Radius filtering: " << m_RadiusFiltering << std::endl;
  os << indent << "Radius estimation: " << m_RadiusEstimation << std::endl;
  os << indent << "Iterations: " << m_Iterations << std::endl;
  os << indent << "Minimum number of used voxels filtering: " << m_MinimumNumberOfUsedVoxelsFiltering << std::endl;
  os << indent << "Use absolute value: " << m_UseAbsoluteValue << std::endl;
  os << indent << "Keep value: " << m_KeepValue << std::endl;

}

} // end namespace itk

#endif
