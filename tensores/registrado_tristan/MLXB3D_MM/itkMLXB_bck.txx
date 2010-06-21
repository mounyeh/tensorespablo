/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMLXB.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkMLXB_txx
#define _itkMLXB_txx
#include "itkMLXB.h"
#include "itkImageFileWriter.h"

namespace itk
{

template < class TInputImage1, class TInputImage2 >
MLXB< TInputImage1, TInputImage2 >::MLXB()
{
	this->SetNumberOfRequiredInputs( 2 );
	this->SetNumberOfRequiredOutputs( 1 );
	
	m_NLevels = 5;
	m_NIter   = 5;
	m_BlockSize.Fill(2);
	m_SearchSize.Fill(1);
	m_Sample.Fill(6);
	
	m_Sigma.SetSize( TInputImage1::ImageDimension );
	m_Sigma.Fill( 16.0 );
	
	m_Metric  = 7;
	m_Bins    = 20;
	
	m_GaussianRadius.Fill( 4 );
	
	m_Output = NULL;
}

/** Start registration */
template < class TInputImage1, class TInputImage2 >
void MLXB< TInputImage1, TInputImage2 >::Start( void )
{
	//----------------------------------------------------------------------------------------------------
	// Cast input images to internal image type:
	FixedCastFilterPointer  fixedCaster   = FixedCastFilterType::New();
	MovingCastFilterPointer movingCaster  = MovingCastFilterType::New();
	fixedCaster->SetInput(  this->GetInput1() );
	movingCaster->SetInput( this->GetInput2() );
	
	// Subsample casted inputs:
	PyramidPointer          fixedPyramid  = PyramidType::New();
	PyramidPointer          movingPyramid = PyramidType::New();
	fixedPyramid->SetInput(  fixedCaster->GetOutput()  );
	movingPyramid->SetInput( movingCaster->GetOutput() );
	
	// Set number of levels in both pyramids:
	fixedPyramid->SetNumberOfLevels(  m_NLevels );
	movingPyramid->SetNumberOfLevels( m_NLevels );
	
	// Update downsampling scheme:
	fixedPyramid->Update();
	movingPyramid->Update();
	//----------------------------------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------------------
	// Create block-matching filter:
	CCFilterPointer         ccFilter;
	// Create bayesian regularization filters:
	RegFilterPointer        regFilter1;
	RegFilterPointer        regFilter2;
	// Create MAP criterion filter:
	MAPFilterPointer        mapFilter;
	// Create transforms:
	TransformPointer        transform1     = TransformType::New();
	TransformPointer        transform2     = TransformType::New();
	
	// Set parameters of the transforms
	transform1->SetGridRegion(    movingCaster->GetOutput()->GetOrigin(), 
		                          movingCaster->GetOutput()->GetSpacing(),
		                          movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
	
	transform1->SetNumberOfPoints( 32, 32, 32 );
	typename TransformType::ParametersType par1;
	par1.SetSize( transform1->GetNumberOfParameters() );
	transform1->SetParameters( par1 );
	transform1->SetIdentity();	
	
	transform2->SetGridRegion(    movingCaster->GetOutput()->GetOrigin(), 
		                          movingCaster->GetOutput()->GetSpacing(),
		                          movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transform2->SetNumberOfPoints( 1, 1, 1 );
	typename TransformType::ParametersType par2;
	par2.SetSize( transform2->GetNumberOfParameters() );
	transform2->SetParameters( par2 );
	transform2->SetIdentity();

	// Create the internal resampling filter for deformation of image 2
	ResampleFilterPointer   resampleFilter;
	InterpolatorPointer     interpolator;
	//----------------------------------------------------------------------------------------------------
	
	for( unsigned int l=0; l<m_NLevels; l++ ){
		this->InvokeEvent( IterationEvent() );
		//------------------------------------------------------------------------------------------------
		// Deform downsampled moving image
		resampleFilter = ResampleFilterType::New();
		interpolator   = InterpolatorType::New();
		resampleFilter->SetInterpolator( interpolator );
		resampleFilter->SetSize(          fixedPyramid->GetOutput( l )->GetLargestPossibleRegion().GetSize() );
		resampleFilter->SetOutputOrigin(  fixedPyramid->GetOutput( l )->GetOrigin()                          );
		resampleFilter->SetOutputSpacing( fixedPyramid->GetOutput( l )->GetSpacing()                         );
		resampleFilter->SetTransform(     transform1                                                         );
		resampleFilter->SetInput(         movingPyramid->GetOutput( l )                                      );
		resampleFilter->Update();
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		// Block-matching for this level of resolution:
		ccFilter       = CCFilterType::New();
		// Set parameters:
		ccFilter->SetMetric(     m_Metric                       );
		ccFilter->SetBins(       m_Bins                         );
		ccFilter->SetSample(     m_Sample                       );
		ccFilter->SetBlockSize(  m_BlockSize                    );
		ccFilter->SetSearchSize( m_SearchSize                   );
		ccFilter->SetInput1(     fixedPyramid->GetOutput( l )   );
		ccFilter->SetInput2(     resampleFilter->GetOutput( )   );
		ccFilter->Update();
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		regFilter1     = RegFilterType::New();
		regFilter2     = RegFilterType::New();
		// Set parameters and initialize:
		regFilter1->SetSearchSize( ccFilter->GetSearchSize() );
		regFilter2->SetSearchSize( ccFilter->GetSearchSize() );
		regFilter1->SetSigma( m_Sigma );
		regFilter2->SetSigma( m_Sigma );
		regFilter1->Initialize();
		regFilter2->Initialize();
		// Bayesian regularization:
		for( int k=0; k<m_NIter; k++ ){
			if( k==0 ){
				regFilter1->SetInput(  ccFilter->GetOutput()   );
				regFilter1->Update();
			}
			else if( k%2 ){
				regFilter2->SetInput(  regFilter1->GetOutput() );
				regFilter1->Update();
			}
			else{
				regFilter1->SetInput(  regFilter1->GetOutput() );
				regFilter1->Update();
			}
		}
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		// MAP criterion:
		mapFilter      = MAPFilterType::New();
		if( m_NIter%2 )
			mapFilter->SetInput( regFilter1->GetOutput() );
		else
			mapFilter->SetInput( regFilter2->GetOutput() );
		mapFilter->Initialize( fixedPyramid->GetOutput( l )->GetSpacing(), ccFilter->GetSearchSize() );
		mapFilter->Update();		
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		SmoothFilterPointer sField = SmoothFilterType::New();
		sField->SetSigma( m_GaussianRadius );
		sField->SetInput( mapFilter->GetOutput() );
		sField->Update();
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		// Interpolate displacement field:
		transform2->SetNumberOfPoints( 1, 1, 1 );
		transform2->SetIdentity();
		
		unsigned int RL = 6;
		for( int k=0; k<RL-1; k++ ){
			transform2->InterpolateImage( sField->GetOutput() );
			transform2->RefineGrid();
		}
		transform2->InterpolateImage(     sField->GetOutput() );
		
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		// Update global deformation field:
		par2.SetSize( transform2->GetNumberOfParameters() );
		par2 = transform2->GetParameters();
		par1 = transform1->GetParameters();
		for( unsigned long p=0; p<transform2->GetNumberOfParameters(); p++ )
			par1[p] += par2[p];
		transform1->SetParameters( par1 );
		//------------------------------------------------------------------------------------------------
		
	}

	//----------------------------------------------------------------------------------------------------
	// Final deformation of the second input:
	ResamplePointer          finalResample     = ResampleType::New();
	FinalInterpolatorPointer finalInterpolator = FinalInterpolatorType::New();
	finalResample->SetInterpolator( finalInterpolator );
	finalResample->SetSize(          this->GetInput2()->GetLargestPossibleRegion().GetSize() );
	finalResample->SetOutputOrigin(  this->GetInput2()->GetOrigin()                          );
	finalResample->SetOutputSpacing( this->GetInput2()->GetSpacing()                         );
	finalResample->SetTransform(     transform1                                              );
	finalResample->SetInput(         this->GetInput2()                                       );
	finalResample->Update();
	m_Output = finalResample->GetOutput();
	//----------------------------------------------------------------------------------------------------

	return;
}

/** Standard "PrintSelf" method */
template < class TInputImage1, class TInputImage2 >
void MLXB< TInputImage1, TInputImage2 >::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << "m_NLevels:           " << m_NLevels    << std::endl;
	os << "m_BlockSize:         " << m_BlockSize  << std::endl;
	os << "m_BlockSize:         " << m_BlockSize  << std::endl;
	os << "m_SearchSize:        " << m_SearchSize << std::endl;
	os << "m_Sample:            " << m_Sample     << std::endl;
	os << "m_Sigma:             " << m_Sigma      << std::endl;
}


} // end namespace itk

#endif
