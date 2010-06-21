/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsStepFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkDemonsStepFilter_txx
#define _itkDemonsStepFilter_txx
#include "itkDemonsStepFilter.h"
#include "itkImageFileWriter.h"

namespace itk
{

template < const unsigned int NDimension >
DemonsStepFilter< NDimension >::DemonsStepFilter()
{
	this->SetNumberOfRequiredInputs(  4  );
	this->SetNumberOfRequiredOutputs( 1  );
	
	m_Lambda             = 1.0;

	m_UseEReg            = true;
	m_UseFReg            = false;
	
	m_SigmaElastic.Fill(   2   );
	m_SigmaFluid.Fill(     2   );
	
	m_SigmaStats         = 4.1410;
	m_SigmaGradient      = 0.577;
	
	m_Field              = NULL;
	
	//----------------------------------------------------------------------------------
	m_LCC                               = ComputeLCCType::New();
	InterpolatePointer      interp      = InterpolateType::New();
	m_Warp                              = WarpFilterType::New();
	m_Warp->SetInterpolator( interp );
	m_LStats                            = LocalStatsFilterType::New();
	m_LJ2Stats                          = LocalStatsFilterType::New();
	m_Square                            = SquareType::New();
	m_LIJStats                          = LocalStatsFilterType::New();
	m_Cross                             = CrossProductType::New();
	m_GradientCalculator                = GradientCalculatorType::New();
	m_Correction                        = ComputeCorrectionType::New();
	m_Add                               = AddFilterType::New();
	m_FluidReg                          = SmoothingFilterType::New();
	m_ElasticReg                        = SmoothingFilterType::New();
	//----------------------------------------------------------------------------------
}

/** Perform a single registration step */
template < const unsigned int NDimension >
void DemonsStepFilter< NDimension >::GenerateData( void )
{
	// Check inputs:
	//     - Input 1: The fixed image I
	//     - Input 2: The smoothed value of I, mean( I )
	//     - Input 3: The smoothed value of I^2, mean( I^2 );
	//     - Input 4: The moving image J
	
	m_LCC->SetInput1( this->GetInput1() );
	m_LCC->SetInput3( this->GetInput2() );
	m_LCC->SetInput5( this->GetInput3() );
	
	// Deform moving image by the current value of the warping filter:
	
	m_Warp->SetInput(            this->GetInput4()               );
	m_Warp->SetOutputSpacing(    this->GetInput1()->GetSpacing() );
	m_Warp->SetOutputOrigin(     this->GetInput1()->GetOrigin()  );
	m_Warp->SetDeformationField( m_Field                         );
	m_Warp->Modified();
	
	m_LCC->SetInput2( m_Warp->GetOutput() );

	m_LStats->SetSigma( m_SigmaStats );
	m_LStats->SetInput( m_Warp->GetOutput()    );
	m_LStats->Modified();
	
	m_LCC->SetInput4(   m_LStats->GetOutput()  );	
	
	// Compute the quadratic of J:
	m_Square->SetInput(    m_Warp->GetOutput()   );
	m_Square->Modified();
	m_LJ2Stats->SetSigma( m_SigmaStats );
	m_LJ2Stats->SetInput(  m_Square->GetOutput() );
	m_LJ2Stats->Modified();
	
	m_LCC->SetInput6( m_LJ2Stats->GetOutput() );
	
	// Compute the cross product IJ:
	m_Cross->SetInput1(    this->GetInput1( )    );
	m_Cross->SetInput2(    m_Warp->GetOutput( )  );
	m_Cross->Modified();
	m_LIJStats->SetSigma( m_SigmaStats );
	m_LIJStats->SetInput(  m_Cross->GetOutput()  );
	m_LIJStats->Modified();
	m_LCC->SetInput7( m_LIJStats->GetOutput() );
	m_LCC->Modified();
	
	// Compute the gradient of the moving image:
	m_GradientCalculator->SetSigma(  m_SigmaGradient*(this->GetInput4()->GetSpacing()[0])  );
	m_GradientCalculator->SetInput(  m_Warp->GetOutput()  );
	m_GradientCalculator->Modified();
	m_GradientCalculator->Update();

	// Compute the correction to the deformation field:
	m_Correction->SetLambda( m_Lambda );
	m_Correction->SetInput1( m_GradientCalculator->GetOutput() );
	m_Correction->SetInput2( m_LCC->GetOutput()                );
	m_Correction->Modified();
	
	// Correct the current deformation field with the (smoothed) velocity field:
	if( m_UseFReg ){
		m_FluidReg->SetSigma(  m_SigmaFluid              );
		m_FluidReg->SetInput(  m_Correction->GetOutput() );
		m_FluidReg->Modified();
		m_Add->SetInput1(      m_FluidReg->GetOutput()   );
	}
	else
		m_Add->SetInput1(   m_Correction->GetOutput()   );
	
	m_Add->SetInput2(   m_Field   );
	m_Add->Modified();
	
	// Execute the whole pipeline:
	if( m_UseEReg ){
		m_ElasticReg->SetSigma(  m_SigmaElastic    );
		m_ElasticReg->SetInput(  m_Add->GetOutput()  );
		m_ElasticReg->Modified();
		m_ElasticReg->Update();
		this->GraftOutput( m_ElasticReg->GetOutput() );
	}
	else{
		m_Add->Update();
		this->GraftOutput( m_Add->GetOutput() );
	}
	return;
}

/** Standard "PrintSelf" method */
template < const unsigned int NDimension >
void DemonsStepFilter< NDimension >::PrintSelf( std::ostream& os, Indent indent ) const
{
	Superclass::PrintSelf( os, indent );
	os << indent << "m_Lambda:        "   <<  m_Lambda         <<  std::endl;
	os << indent << "m_SigmaElastic:  "   <<  m_SigmaElastic   <<  std::endl;
	os << indent << "m_SigmaFluid:    "   <<  m_SigmaFluid     <<  std::endl;
	os << indent << "m_SigmaStats:    "   <<  m_SigmaStats     <<  std::endl;
	os << indent << "m_SigmaGradient: "   <<  m_SigmaGradient  <<  std::endl;
}


} // end namespace itk

#endif
