/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMLD.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkMLD_txx
#define _itkMLD_txx
#include "itkMLD.h"
//===================================================================================================
//===================================================================================================
//===================================================================================================
// ESTE BLOQUE ES SOLO PARA VALIDACIÓN!!!
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
//===================================================================================================
//===================================================================================================
//===================================================================================================

namespace itk
{

template < class TInputImage1, class TInputImage2 >
MLD< TInputImage1, TInputImage2 >::MLD()
{
	this->SetNumberOfRequiredInputs(  2  );
	this->SetNumberOfRequiredOutputs( 1  );
	
	m_NLevels            = 5;
	m_Steps              = 10;
	
	m_Lambda             = 1.0;
	m_Tau                = exp(  -0.693147/((double)m_Steps)  );

	m_SigmaElastic.Fill(  2  );
	m_SigmaFluid.Fill(    9  );
	
	m_SigmaStats         = 4.1410;
	m_SigmaGradient      = 0.577;

    m_UseEReg            = true;
    m_UseFReg            = true;

	m_Output = NULL;
	m_DeformationField = NULL;
}

/** Start registration */
template < class TInputImage1, class TInputImage2 >
void MLD< TInputImage1, TInputImage2 >::Start( void )
{
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	// Not-tested code
	if(   static_cast<DeformationFieldType*>( m_DeformationField )   !=   NULL   )
		m_NLevels = 1;
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------------------
	// Check if the number of levels is excesive or if the images to register has different sizes:
	unsigned long minSize = this->GetInput1()->GetLargestPossibleRegion().GetSize()[0];
	for( unsigned int dim=0; dim<TInputImage1::ImageDimension; dim++ ){
		if( this->GetInput1()->GetLargestPossibleRegion().GetSize()[dim] < minSize )
			minSize = this->GetInput1()->GetLargestPossibleRegion().GetSize()[dim];
		if(   this->GetInput1()->GetRequestedRegion().GetSize()[dim] != this->GetInput2()->GetRequestedRegion().GetSize()[dim]   )
			itkExceptionMacro( << "Requested regions of fixed and moving image should be the same size" );
	}
	double minResolution = (double)minSize;
	for( unsigned int r=1; r<m_NLevels; r++ )
		minResolution /= 2.0;
	while( minResolution<10.0 ){
		m_NLevels--;
		minResolution *= 2.0;
	}
	//----------------------------------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------------------
	// Cast input images to internal image type:
	FixedCastFilterPointer  fixedCaster   = FixedCastFilterType::New();
	MovingCastFilterPointer movingCaster  = MovingCastFilterType::New();
	fixedCaster->SetInput(  this->GetInput1() );
	movingCaster->SetInput( this->GetInput2() );
	
	// Subsample casted inputs:
	PyramidPointer          fixedPyramid  = PyramidType::New();
	PyramidPointer          movingPyramid = PyramidType::New();
	fixedPyramid->SetInput(    fixedCaster->GetOutput()    );
	movingPyramid->SetInput(   movingCaster->GetOutput()   );
	
	// Set number of levels in both pyramids:
	fixedPyramid->SetNumberOfLevels(  m_NLevels );
	movingPyramid->SetNumberOfLevels( m_NLevels );
	
	// Update the downsampling scheme:
	fixedPyramid->Update();
	movingPyramid->Update();
	//----------------------------------------------------------------------------------------------------

	//----------------------------------------------------------------------------------------------------
	DeformationFieldTypePointer field;
	
	// For each resolution level:
	for( unsigned int level=0; level<m_NLevels; level++ ){
		this->InvokeEvent( IterationEvent() );
		// Creation of the Demons step object:
		DemonsPointer demons;
		// Get the square of the fixed image:
		SquarePointer square = SquareType::New();
		square->SetInput(   fixedPyramid->GetOutput( level )   );
		// Create local means filters:
		LocalStatsFilterPointer localStats  = LocalStatsFilterType::New();
		LocalStatsFilterPointer qlocalStats = LocalStatsFilterType::New();
		localStats->SetSigma(   m_SigmaStats   );
		qlocalStats->SetSigma(  m_SigmaStats   );
		localStats->SetInput(   fixedPyramid->GetOutput( level )   );
		qlocalStats->SetInput(  square->GetOutput( )               );
		localStats->Update();
		qlocalStats->Update();
		// Value of the Levenberg-Marquardt parameter in each iteration:
		double lambda = m_Lambda;
		//-------------------------------------------------------------------------------------------------
		if( level==0 ){
			if(   static_cast<DeformationFieldType*>( m_DeformationField )   ==   NULL   ){// Not-tested code
				// If this is the first level, we have to create the field and the temporal field:
				typename DeformationFieldType::PixelType initialValue( 0.0 );
				field = DeformationFieldType::New();
				field->SetRegions(      fixedPyramid->GetOutput( 0 )->GetLargestPossibleRegion()  );
				field->Allocate();
				field->FillBuffer(      initialValue  );
			}
			else{// Not-tested code
				field = m_DeformationField;// Not-tested code
			}
		}
		//-------------------------------------------------------------------------------------------------
		// Iterations:
		demons = DemonsType::New();
		for( unsigned int step=0; step<m_Steps; step++ ){
			// Set the currently available inputs:
			demons->SetUseFluidRegularization(   this->GetUseFluidRegularization()   );
			demons->SetUseElasticRegularization( this->GetUseElasticRegularization() );
			demons->SetInput1(   fixedPyramid->GetOutput( level )    );
			demons->SetInput4(   movingPyramid->GetOutput( level )   );
			// Get a smoothed version of the fixed image and set the 3rd input to the demons step:
			demons->SetInput2(      localStats->GetOutput()      );
			// Get a smoothed version of the squared fixed image and set the 4th input to the demons step:
			demons->SetInput3(      qlocalStats->GetOutput()     );
			// Parameter tunning:
			demons->SetSigmaElastic(  m_SigmaElastic  );
			demons->SetSigmaFluid(    m_SigmaFluid    );
			demons->SetSigmaStats(    m_SigmaStats    );
			demons->SetSigmaGradient( m_SigmaGradient );
			// Current value of lambda:
			demons->SetLambda( lambda );
			// Set the appropriate 5th input:
			demons->SetDeformationField( field );
			// Single correction step
			demons->Update();
			// Get the new field:
			field = demons->GetOutput();
			// Correction of lambda:
			lambda *= m_Tau;
		}
				
		// If this is not the last resolution level, we must obtain a new deformation field from the
		// previous decimated version
		/** Bloque comprobado (funciona bien, y field queda bien escalado */
		if( level != m_NLevels-1 ){
			// Create a new upsampling filter:
			ResampleFilterPointer upsample = ResampleFilterType::New();
			// Set the output spacing and origin:
			upsample->SetOutputSpacing( movingPyramid->GetOutput( level+1 )->GetSpacing() );
			upsample->SetOutputOrigin(  movingPyramid->GetOutput( level+1 )->GetOrigin()  );
			// Set the new size for the deformation field:
			upsample->SetSize( movingPyramid->GetOutput( level+1 )->GetLargestPossibleRegion().GetSize() );
			// Set the appropriate input:
			upsample->SetInput( field );
			// Update the filter:
			upsample->UpdateLargestPossibleRegion();
			// Get the upsampled field:
			field = upsample->GetOutput();
			// Destroy the filter:
			upsample = NULL;
		}
	}
	//----------------------------------------------------------------------------------------------------
	// Final deformation of the second input:
	WarpFilterPointer   warp   = WarpFilterType::New();
	InterpolatorPointer interp = InterpolatorType::New();
	// Set the parameters of the warping filter:
	warp->SetInput(              m_Input2                 );
	warp->SetInterpolator(       interp                   );
	warp->SetOutputSpacing(      m_Input1->GetSpacing( )  );
	warp->SetOutputOrigin(       m_Input1->GetOrigin(  )  );
	warp->SetDeformationField(   field                    );
	// Update the filter:
	warp->Update();
	m_Output = warp->GetOutput();
	// Set the new deformation field:
	m_DeformationField = field;
	field = NULL;
	//----------------------------------------------------------------------------------------------------
	
	
	return;
}

/** Standard "PrintSelf" method */
template < class TInputImage1, class TInputImage2 >
void MLD< TInputImage1, TInputImage2 >::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "m_NLevels:       "   <<  m_NLevels        <<  std::endl;
	os << indent << "m_Steps:         "   <<  m_Steps          <<  std::endl;
	os << indent << "m_Lambda:        "   <<  m_Lambda         <<  std::endl;
	os << indent << "m_Tau:           "   <<  m_Tau            <<  std::endl;
	os << indent << "m_SigmaElastic:  "   <<  m_SigmaElastic   <<  std::endl;
	os << indent << "m_SigmaFluid:    "   <<  m_SigmaFluid     <<  std::endl;
	os << indent << "m_SigmaStats:    "   <<  m_SigmaStats     <<  std::endl;
	os << indent << "m_SigmaGradient: "   <<  m_SigmaGradient  <<  std::endl;
}


} // end namespace itk

#endif
