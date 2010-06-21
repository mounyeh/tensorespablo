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
	m_BlockSize.Fill(5);
	m_SearchSize.Fill(2);
	m_Sample.Fill(6);
	
	m_Sigma.SetSize( TInputImage1::ImageDimension );
	m_Sigma.Fill( 16.0 );
	
	m_Lambda  = 0.5;
	
	m_GaussianRadius.Fill( 2 );
	
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
	// The number of refinement level in the interpoaltion:
	unsigned int RL = 5;
	unsigned int NP = 1;
	for( int k=0; k<RL; k++ )
		NP*=2;
	//----------------------------------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------------------
	// Create block-matching filters pair:
	CCFilterPointer         ccFilterf;
	CCFilterPointer         ccFilteri;
	// Create bayesian regularization filters:
	RegFilterPointer        regFilter1f;
	RegFilterPointer        regFilter2f;
	RegFilterPointer        regFilter1i;
	RegFilterPointer        regFilter2i;
	// Create MAP criterion filter:
	MAPFilterPointer        mapFilterf;
	MAPFilterPointer        mapFilteri;
	// Create transforms:
	TransformPointer        transform1      = TransformType::New();
	TransformPointer        transformf2     = TransformType::New();
	TransformPointer        transformi2     = TransformType::New();
	
	// Set parameters of the transforms
	transform1->SetGridRegion(    movingCaster->GetOutput()->GetOrigin(), 
		                          movingCaster->GetOutput()->GetSpacing(),
		                          movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transform1->SetNumberOfPoints( NP, NP, NP );
	
	typename TransformType::ParametersType par1;
	par1.SetSize( transform1->GetNumberOfParameters() );
	transform1->SetParameters( par1 );
	transform1->SetIdentity();
	
	transformf2->SetGridRegion(   movingCaster->GetOutput()->GetOrigin(), 
		                          movingCaster->GetOutput()->GetSpacing(),
		                          movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transformf2->SetNumberOfPoints( 1, 1, 1 );
	transformi2->SetGridRegion(   movingCaster->GetOutput()->GetOrigin(), 
		                          movingCaster->GetOutput()->GetSpacing(),
		                          movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transformi2->SetNumberOfPoints( 1, 1, 1 );
	typename TransformType::ParametersType parf2;
	parf2.SetSize( transformf2->GetNumberOfParameters() );
	typename TransformType::ParametersType pari2;
	pari2.SetSize( transformi2->GetNumberOfParameters() );
	transformf2->SetParameters( parf2 );
	transformf2->SetIdentity();
	transformi2->SetParameters( pari2 );
	transformi2->SetIdentity();

	// Create the internal resampling filter for deformations
	ResampleFilterPointer   resampleFilter;
	InterpolatorPointer     interpolator;
	//----------------------------------------------------------------------------------------------------


	for( unsigned int l=0; l<m_NLevels; l++ ){
		this->InvokeEvent( IterationEvent() );
		//------------------------------------------------------------------------------------------------
		// Deform downsampled moving image
		resampleFilter  = ResampleFilterType::New();
		interpolator    = InterpolatorType::New();
		resampleFilter->SetInterpolator( interpolator );
		resampleFilter->SetSize(          movingPyramid->GetOutput( l )->GetLargestPossibleRegion().GetSize() );
		resampleFilter->SetOutputOrigin(  movingPyramid->GetOutput( l )->GetOrigin()                          );
		resampleFilter->SetOutputSpacing( movingPyramid->GetOutput( l )->GetSpacing()                         );
		resampleFilter->SetTransform(     transform1                                                          );
		resampleFilter->SetInput(         movingPyramid->GetOutput( l )                                       );
		resampleFilter->Update();
		//------------------------------------------------------------------------------------------------

		//------------------------------------------------------------------------------------------------
		// Block-matching for this level of resolution:
		// --------------- The forward block-matching:
		ccFilterf = CCFilterType::New();
		// Set parameters:
		ccFilterf->SetSample(     m_Sample                       );
		ccFilterf->SetBlockSize(  m_BlockSize                    );
		ccFilterf->SetSearchSize( m_SearchSize                   );
		ccFilterf->SetInput1(     fixedPyramid->GetOutput( l )   );
		ccFilterf->SetInput2(     resampleFilter->GetOutput( )   );
		ccFilterf->Update();
		// --------------- The inverse block-matching:
		ccFilteri = CCFilterType::New();
		// Set parameters:
		ccFilteri->SetSample(     m_Sample                       );
		ccFilteri->SetBlockSize(  m_BlockSize                    );
		ccFilteri->SetSearchSize( m_SearchSize                   );
		ccFilteri->SetInput1(     resampleFilter->GetOutput( )   );
		ccFilteri->SetInput2(     fixedPyramid->GetOutput( l )   );
		ccFilteri->Update();
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		// ----------------- Regularization of the forward deformation field:
		regFilter1f = RegFilterType::New();
		regFilter2f = RegFilterType::New();
		// Set parameters and initialize:
		regFilter1f->SetSearchSize( ccFilterf->GetSearchSize() );
		regFilter2f->SetSearchSize( ccFilterf->GetSearchSize() );
		regFilter1f->SetSigma( m_Sigma );
		regFilter2f->SetSigma( m_Sigma );
		regFilter1f->Initialize();
		regFilter2f->Initialize();
		// Bayesian regularization:
		for( int k=0; k<m_NIter; k++ ){
			if( k==0 ){
				regFilter1f->SetInput(  ccFilterf->GetOutput()   );
				regFilter1f->Update();
			}
			else if( k%2 ){
				regFilter2f->SetInput(  regFilter1f->GetOutput() );
				regFilter2f->Update();
			}
			else{
				regFilter1f->SetInput(  regFilter2f->GetOutput() );
				regFilter1f->Update();
			}
		}
		// ----------------- Regularization of the inverse deformation field:
		regFilter1i = RegFilterType::New();
		regFilter2i = RegFilterType::New();
		// Set parameters and initialize:
		regFilter1i->SetSearchSize( ccFilteri->GetSearchSize() );
		regFilter2i->SetSearchSize( ccFilteri->GetSearchSize() );
		regFilter1i->SetSigma( m_Sigma );
		regFilter2i->SetSigma( m_Sigma );
		regFilter1i->Initialize();
		regFilter2i->Initialize();
		// Bayesian regularization:
		for( int k=0; k<m_NIter; k++ ){
			if( k==0 ){
				regFilter1i->SetInput(  ccFilteri->GetOutput()   );
				regFilter1i->Update();
			}
			else if( k%2 ){
				regFilter2i->SetInput(  regFilter1i->GetOutput() );
				regFilter2i->Update();
			}
			else{
				regFilter1i->SetInput(  regFilter2i->GetOutput() );
				regFilter1i->Update();
			}
		}
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		// ------------------ MAP criterion for the forward transform:
		mapFilterf = MAPFilterType::New();
		if( m_NIter==0 )
			mapFilterf->SetInput( ccFilterf->GetOutput() );
		else if( m_NIter%2 )
			mapFilterf->SetInput( regFilter1f->GetOutput() );
		else
			mapFilterf->SetInput( regFilter2f->GetOutput() );
		mapFilterf->Initialize( fixedPyramid->GetOutput( l )->GetSpacing(), ccFilterf->GetSearchSize() );
		mapFilterf->Update();
		// ------------------ MAP criterion for the inverse transform:
		mapFilteri = MAPFilterType::New();
		if( m_NIter==0 )
			mapFilteri->SetInput( ccFilteri->GetOutput() );
		else if( m_NIter%2 )
			mapFilteri->SetInput( regFilter1i->GetOutput() );
		else
			mapFilteri->SetInput( regFilter2i->GetOutput() );
		mapFilteri->Initialize( movingPyramid->GetOutput( l )->GetSpacing(), ccFilteri->GetSearchSize() );
		mapFilteri->Update();		
		//------------------------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------------------------
		// ----------------- Gaussian smoothing of the forward transform:
		SmoothFilterPointer sFieldf = SmoothFilterType::New();
		sFieldf->SetSigma( m_GaussianRadius );
		sFieldf->SetInput( mapFilterf->GetOutput() );
		sFieldf->Update();
		// ----------------- Gaussian smoothing of the inverse transform:
		SmoothFilterPointer sFieldi = SmoothFilterType::New();
		sFieldi->SetSigma( m_GaussianRadius );
		sFieldi->SetInput( mapFilteri->GetOutput() );
		sFieldi->Update();
		//------------------------------------------------------------------------------------------------

		//------------------------------------------------------------------------------------------------
		// ----------------------------- Interpolate forward displacement field:
		transformf2->SetNumberOfPoints( 1, 1, 1 );
		transformf2->SetGridRegion(   movingCaster->GetOutput()->GetOrigin(), 
					      movingCaster->GetOutput()->GetSpacing(),
		                              movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
		transformf2->SetIdentity();
		
		for( int k=0; k<RL; k++ ){
			transformf2->InterpolateImage( sFieldf->GetOutput() );
			transformf2->RefineGrid();
		}
		transformf2->InterpolateImage(     sFieldf->GetOutput() );
		// ----------------------------- Interpolate inverse displacement field:
		transformi2->SetNumberOfPoints( 1, 1, 1 );
		transformi2->SetGridRegion(   movingCaster->GetOutput()->GetOrigin(), 
					      movingCaster->GetOutput()->GetSpacing(),
		                              movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
		transformi2->SetIdentity();
		
		for( int k=0; k<RL; k++ ){
			transformi2->InterpolateImage( sFieldi->GetOutput() );
			transformi2->RefineGrid();
		}
		transformi2->InterpolateImage(     sFieldi->GetOutput() );
		//------------------------------------------------------------------------------------------------
		
		// -------------------- Update global deformation field for forward transform:
		parf2.SetSize( transformf2->GetNumberOfParameters() );
		pari2.SetSize( transformi2->GetNumberOfParameters() );
		parf2 = transformf2->GetParameters();
		pari2 = transformi2->GetInverseParameters();
		par1  = transform1->GetParameters();
		for( unsigned long p=0; p<transform1->GetNumberOfParameters(); p++ )
			par1[p] += m_Lambda*pari2[p] + (1.0-m_Lambda)*parf2[p];
		transform1->SetParameters( par1 );
		//------------------------------------------------------------------------------------------------
		
		
		
	}

	//----------------------------------------------------------------------------------------------------
	//std::cout << transform1 << std::endl;
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
	os << "m_SearchSize:        " << m_SearchSize << std::endl;
	os << "m_Sample:            " << m_Sample     << std::endl;
	os << "m_Sigma:             " << m_Sigma      << std::endl;
}


} // end namespace itk

#endif
