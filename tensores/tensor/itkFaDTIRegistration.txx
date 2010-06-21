#ifndef _itkFaDTIRegistration_txx
#define _itkFaDTIRegistration_txx

#include "itkFaDTIRegistration.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

template < class TInputImage1, class TInputImage2, class TLabeledImage >
FaDTIRegistration< TInputImage1, TInputImage2, TLabeledImage >::FaDTIRegistration()
{
	this->SetNumberOfRequiredInputs( 3 );
	this->SetNumberOfRequiredOutputs( 1 );
	
	m_NLevels = 4;
	m_NIter   = 5;
	
	m_Sample.Fill(2);
	m_BlockSize.Fill(1);
	m_SearchSize.Fill(1);

	m_Lambda  = 0.5;
	
	m_GaussianRadius.Fill( 2 );

	m_Sigma.SetSize( TInputImage1::ImageDimension );
	m_Sigma.Fill( 16.0 );

	m_Output = NULL;
}

/** Start registration */
template < class TInputImage1, class TInputImage2, class TLabeledImage >
void FaDTIRegistration< TInputImage1, TInputImage2, TLabeledImage >::Start( void )
{

	
	ComputeScalarsPointer fixedFAFilter		=	ComputeScalarsType::New();
	ComputeScalarsPointer movingFAFilter	=	ComputeScalarsType::New();

	FixedImagePointer          fixedT		=	FixedImageType::New();
	FixedImagePointer   	   movingT		=	FixedImageType::New();
	
	LabelImagePointer	labels  = LabelImageType::New();
	
	fixedT=this->GetFixed();
	
	labels=this->GetLabels();
	labels->SetOrigin(fixedT->GetOrigin());
	labels->SetSpacing(fixedT->GetSpacing());
	
	movingT=this->GetMoving();
	movingT->SetOrigin(fixedT->GetOrigin());
	movingT->SetSpacing(fixedT->GetSpacing());
	
	
	//fixedFAFilter->SetInput( this->GetFixed() );
	//movingFAFilter->SetInput( this->GetMoving() );
	fixedFAFilter->SetInput( fixedT );
	movingFAFilter->SetInput( movingT );


	fixedFAFilter->SetComputeFA();
	movingFAFilter->SetComputeFA();

	try{
		fixedFAFilter->Update();
		movingFAFilter->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}

	//mask	--> true if we are in the background of both image
	//		--> false elsewhere
	MaskImageType::Pointer mask	=MaskImageType::New();
	

std::cout<<fixedFAFilter->GetOutput()->GetOrigin()<<std::endl;
std::cout<<fixedFAFilter->GetOutput()->GetSpacing()<<std::endl;
std::cout<<movingFAFilter->GetOutput()->GetOrigin()<<std::endl;
std::cout<<movingFAFilter->GetOutput()->GetSpacing()<<std::endl;
std::cout<<this->GetLabels()->GetOrigin()<<std::endl;
std::cout<<this->GetLabels()->GetSpacing()<<std::endl;
				
	PyramidPointer          fixedPyramid  = PyramidType::New();
	PyramidPointer          movingPyramid = PyramidType::New();
	fixedPyramid->SetInput(  fixedFAFilter->GetOutput()  );
	movingPyramid->SetInput( movingFAFilter->GetOutput() );
	
	//MaskPyramidPointer         maskPyramid = MaskPyramidType::New();
	//maskPyramid->SetInput( mask  );
	
	// Set number of levels in both pyramids:
	fixedPyramid->SetNumberOfLevels(  m_NLevels );
	movingPyramid->SetNumberOfLevels( m_NLevels );
	
	//maskPyramid->SetNumberOfLevels( m_NLevels );
	
	unsigned int shrink[3]={1,1,1};
	SizeType size=(this->GetFixed()->GetRequestedRegion().GetSize());
	
	for(unsigned i=0; i<(m_NLevels-1); i++){
	  if(shrink[0]<(size[0]/8)){
	    shrink[0]= shrink[0]*2;
	  }
	  if(shrink[1]<(size[1]/8)){
	    shrink[1]= shrink[1]*2;
	  }
	  if(shrink[2]<(size[2]/8)){
	    shrink[2]= shrink[2]*2;
	  }
	}
	  
	fixedPyramid->SetStartingShrinkFactors(shrink);
	movingPyramid->SetStartingShrinkFactors(shrink);
//	maskPyramid->SetStartingShrinkFactors(shrink);
	
	// Update downsampling scheme:
	fixedPyramid->Update();
	movingPyramid->Update();
//	maskPyramid->Update();
	
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
	//CCFilterPointer         ccFilteri;
	// Create bayesian regularization filters:
	RegFilterPointer        regFilter1f;
	RegFilterPointer        regFilter2f;
	//RegFilterPointer        regFilter1i;
	//RegFilterPointer        regFilter2i;
	// Create MAP criterion filter:
	MAPFilterPointer        mapFilterf;
	//MAPFilterPointer        mapFilteri;
	// Create transforms:
	TransformPointer        transform1      = TransformType::New();
	TransformPointer        transformf2     = TransformType::New();
	//TransformPointer        transformi2     = TransformType::New();

	// Set parameters of the transforms
	transform1->SetGridRegion(    movingFAFilter->GetOutput()->GetOrigin(), 
		                          movingFAFilter->GetOutput()->GetSpacing(),
		                          movingFAFilter->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transform1->SetNumberOfPoints( NP, NP, NP );
	
	//typename TransformType::ParametersType par1;
	m_Parameters.SetSize( transform1->GetNumberOfParameters() );
	transform1->SetParameters( m_Parameters );
	transform1->SetIdentity();
	
	transformf2->SetGridRegion(   movingFAFilter->GetOutput()->GetOrigin(), 
		                          movingFAFilter->GetOutput()->GetSpacing(),
		                          movingFAFilter->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transformf2->SetNumberOfPoints( 1, 1, 1 );
//	transformi2->SetGridRegion(   movingFAFilter->GetOutput()->GetOrigin(), 
//		                          movingFAFilter->GetOutput()->GetSpacing(),
//		                          movingFAFilter->GetOutput()->GetLargestPossibleRegion().GetSize() );
//	transformi2->SetNumberOfPoints( 1, 1, 1 );
	typename TransformType::ParametersType parf2;
	parf2.SetSize( transformf2->GetNumberOfParameters() );
//	typename TransformType::ParametersType pari2;
//	pari2.SetSize( transformi2->GetNumberOfParameters() );
	transformf2->SetParameters( parf2 );
	transformf2->SetIdentity();
//	transformi2->SetParameters( pari2 );
//	transformi2->SetIdentity();

	// Create the internal resampling filter for deformations
	ResampleFilterPointer   resampleFilter;
	InterpolatorPointer     interpolator;
	//----------------------------------------------------------------------------------------------------

	
	typedef itk::ImageRegionConstIterator<InternalImageType>		ImageIteratorType;
	typedef itk::ImageRegionIterator<MaskImageType>					MaskIteratorType;
	
	for( unsigned int l=0; l<m_NLevels; l++ ){
		this->InvokeEvent( IterationEvent() );
		mask=NULL;
		mask=MaskImageType::New();
		mask->SetRegions(movingPyramid->GetOutput(l)->GetRequestedRegion());
		mask->SetOrigin(movingPyramid->GetOutput(l)->GetOrigin());
		mask->SetSpacing(movingPyramid->GetOutput(l)->GetSpacing());
		mask->Allocate();
		if (l>(m_NLevels-2)){ 
			mask->FillBuffer(true);
			ImageIteratorType it_fix(fixedPyramid->GetOutput(l),fixedPyramid->GetOutput(l)->GetRequestedRegion()); 
			ImageIteratorType it_mov(movingPyramid->GetOutput(l),movingPyramid->GetOutput(l)->GetRequestedRegion()); 
	
			MaskIteratorType it_mask(mask,mask->GetRequestedRegion());	

			for(it_fix.GoToBegin(), it_mov.GoToBegin(), it_mask.GoToBegin(); !it_mov.IsAtEnd(); ++it_fix, ++it_mov, ++it_mask){
				if(it_fix.Get()!=0 || it_mov.Get()!=0){
					it_mask.Set(false);
				}
			}
		}else{
				mask->FillBuffer(false);
		}


		//------------------------------------------------------------------------------------------------
		// Deform downsampled moving image
		resampleFilter  = ResampleFilterType::New();
		interpolator    = InterpolatorType::New();
		 //interpolator->SetMask(mask );
		resampleFilter->SetMask(mask);
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
		ccFilterf->SetMask(mask);
		ccFilterf->SetSample(     m_Sample                       );
		ccFilterf->SetBlockSize(  m_BlockSize                    );
		ccFilterf->SetSearchSize( m_SearchSize                   );
		ccFilterf->SetInput1(     fixedPyramid->GetOutput( l )   );
		ccFilterf->SetInput2(     resampleFilter->GetOutput( )   );
		ccFilterf->Update();
		// --------------- The inverse block-matching:
	//	ccFilteri = CCFilterType::New();
		// Set parameters:
	//	ccFilteri->SetSample(     m_Sample                       );
	//	ccFilteri->SetBlockSize(  m_BlockSize                    );
	//	ccFilteri->SetSearchSize( m_SearchSize                   );
	//	ccFilteri->SetInput1(     resampleFilter->GetOutput( )   );
	//	ccFilteri->SetInput2(     fixedPyramid->GetOutput( l )   );
	//	ccFilteri->Update();
		//------------------------------------------------------------------------------------------------
	
		//------------------------------------------------------------------------------------------------
		
		// ----------------- Regularization of the forward deformation field:
		regFilter1f = RegFilterType::New();
		regFilter2f = RegFilterType::New();
		// Set parameters and initialize:
		regFilter1f->SetSearchSize( ccFilterf->GetSearchSize() );
		regFilter2f->SetSearchSize( ccFilterf->GetSearchSize() );
		regFilter1f->SetMask( mask );
		regFilter2f->SetMask( mask );
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
	/*	regFilter1i = RegFilterType::New();
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
	*/

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
	/*	mapFilteri = MAPFilterType::New();
		if( m_NIter==0 )
			mapFilteri->SetInput( ccFilteri->GetOutput() );
		else if( m_NIter%2 )
			mapFilteri->SetInput( regFilter1i->GetOutput() );
		else
			mapFilteri->SetInput( regFilter2i->GetOutput() );
		mapFilteri->Initialize( movingPyramid->GetOutput( l )->GetSpacing(), ccFilteri->GetSearchSize() );
		mapFilteri->Update();		
		//------------------------------------------------------------------------------------------------
*/
		//------------------------------------------------------------------------------------------------
		// ----------------- Gaussian smoothing of the forward transform:
		SmoothFilterPointer sFieldf = SmoothFilterType::New();
		sFieldf->SetSigma( m_GaussianRadius );
		sFieldf->SetInput( mapFilterf->GetOutput() );
		sFieldf->Update();
		// ----------------- Gaussian smoothing of the inverse transform:
/*		SmoothFilterPointer sFieldi = SmoothFilterType::New();
		sFieldi->SetSigma( m_GaussianRadius );
		sFieldi->SetInput( mapFilteri->GetOutput() );
		sFieldi->Update();
		//------------------------------------------------------------------------------------------------
*/
		//------------------------------------------------------------------------------------------------
		// ----------------------------- Interpolate forward displacement field:
		transformf2->SetNumberOfPoints( 1, 1, 1 );
		transformf2->SetGridRegion(   movingFAFilter->GetOutput()->GetOrigin(), 
									  movingFAFilter->GetOutput()->GetSpacing(),
		                              movingFAFilter->GetOutput()->GetLargestPossibleRegion().GetSize() );
		transformf2->SetIdentity();
		
		for( int k=0; k<RL; k++ ){
			transformf2->InterpolateImage( sFieldf->GetOutput() );
			transformf2->RefineGrid();
		}
		transformf2->InterpolateImage(     sFieldf->GetOutput() );
		// ----------------------------- Interpolate inverse displacement field:
/*		transformi2->SetNumberOfPoints( 1, 1, 1 );
		transformi2->SetGridRegion(   movingFAFilter->GetOutput()->GetOrigin(), 
									  movingFAFilter->GetOutput()->GetSpacing(),
		                              movingFAFilter->GetOutput()->GetLargestPossibleRegion().GetSize() );
		transformi2->SetIdentity();
		
		for( int k=0; k<RL; k++ ){
			transformi2->InterpolateImage( sFieldi->GetOutput() );
			transformi2->RefineGrid();
		}
		transformi2->InterpolateImage(     sFieldi->GetOutput() );
		//------------------------------------------------------------------------------------------------
*/		
		// -------------------- Update global deformation field for forward transform:
		parf2.SetSize( transformf2->GetNumberOfParameters() );
//		pari2.SetSize( transformi2->GetNumberOfParameters() );
		parf2 = transformf2->GetParameters();
//		pari2 = transformi2->GetInverseParameters();
		m_Parameters  = transform1->GetParameters();
		for( unsigned long p=0; p<transform1->GetNumberOfParameters(); p++ )
//			par1[p] += m_Lambda*pari2[p] + (1.0-m_Lambda)*parf2[p];
			m_Parameters[p] += parf2[p];
		
		transform1->SetParameters( m_Parameters );
		//------------------------------------------------------------------------------------------------
	}
		// Final deformation of the second input:

	FinalResamplePointer          finalResample		=	FinalResampleType::New();
	FinalInterpolatorPointer finalInterpolator      =	FinalInterpolatorType::New();

	finalResample->SetInterpolator( finalInterpolator );
	finalResample->SetSize(         labels->GetLargestPossibleRegion().GetSize() );
	finalResample->SetOutputOrigin( labels->GetOrigin()                         );
	finalResample->SetOutputSpacing(labels->GetSpacing()                        );
	finalResample->SetTransform(     transform1                                             );
	finalResample->SetMask(mask);
	finalResample->SetInput(       labels					      				);

	finalResample->Update();
	m_Output=finalResample->GetOutput();
	
	//----------------------------------------------------------------------------------------------------

	return;


}	

/** Standard "PrintSelf" method */
template < class TInputImage1, class TInputImage2, class TLabeledImage >
void FaDTIRegistration< TInputImage1, TInputImage2, TLabeledImage >::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
}


}//end namespace itk

#endif
