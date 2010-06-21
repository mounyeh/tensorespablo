/*=========================================================================

  Program:   FilteringConsole.cxx
  Language:  C++
  Date:      28-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#include "FilteringGUI.h"
#include "FilteringConsole.h"

// itk includes 
#include <itkImageRegionIteratorWithIndex.h>
#include <itkVTKImageExport.h>
// vtk includes 
#include <vtkITKUtility.h>
#include <vtkImageImport.h>
// std includes
#include <vector>

FilteringConsole::FilteringConsole(int X, int Y, int W, int H, const char *label)
:FilteringGUI(X,Y,W,H,label) 
{
	m_tipo_filtroWiener->add("Homogeneous",0,NULL,NULL,0);
	m_tipo_filtroWiener->add("Anisotropic",0,NULL,NULL,0);
	m_tipo_filtroWiener->add("Mixed",0,NULL,NULL,0);
}

FilteringConsole::~FilteringConsole()
{
}

void FilteringConsole::SetModelDataBrowser(Fl_Browser *modeldataBrowserIn)
{
	m_modeldataBrowser = modeldataBrowserIn;
}

void FilteringConsole::SetImageViewer3D(Viewer3D* ImageViewer3DIn)
{
	ImageViewer3D = ImageViewer3DIn;
}

void FilteringConsole::SetVectorModelData(void* VectorModelData)
{
	m_VectorModelData = (VectorOfModelType*)VectorModelData;
}

void FilteringConsole::SetVectorData(void* VectorData)
{
	m_VectorData = (VectorOfDataType*)VectorData;
}

void FilteringConsole::SetImageColorViewer(MyColorViewerType* ImageColorViewerIn)
{
	ImageColorViewer = ImageColorViewerIn;
}

void FilteringConsole::SetSliders(Fl_Value_Slider **sliders)
{
	sliceNumberSlider = sliders;
}

void FilteringConsole::SetProgressSliders(fltk::ProgressBar *progressbar)
{
	progressSlider = progressbar;
}

void FilteringConsole::SetProgressCounter(Fl_Value_Output *progressCounter)
{
	m_ProgressCounter = progressCounter;
}

void FilteringConsole::SetThreads(int threads)
{
	m_threads = threads;
}

void FilteringConsole::SetNImagesCargadas(unsigned int &NImagesCargadas)
{
	m_NImagesCargadas = NImagesCargadas;
}

void FilteringConsole::OnTipoFiltradoChange()
{
	m_tipofiltrado = m_filtrado->value();
}

void FilteringConsole::OnTipoFiltroWienerChange()
{
	m_filtroWiener = m_tipo_filtroWiener->value();
}


void FilteringConsole::OnProgress(itk::Object *object, const itk::EventObject &event)
{
	// Get the value of the progress
	itk::ProcessObject *po = reinterpret_cast<ProcessObject *>(object);
	float progress = po->GetProgress();
	
	//std::cout << "progress "  << progress << std::endl;
	
	float max = 1.0;
	m_ProgressCounter->value(100*progress/max);
	
	// Display the filter's progress
	progressSlider->value(100*progress/max);
	
	// Let the UI refresh
	Fl::check();
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::BilateralFilter(void ) {
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &FilteringConsole::OnProgress );
	
	double DomainSigma[3];
	DomainSigma[0] = BilatDomainSigma0->value();
	DomainSigma[1] = BilatDomainSigma1->value();
	DomainSigma[2] = BilatDomainSigma2->value();
	
	BilateralFilterType::Pointer bilateralFilter = BilateralFilterType::New();
	bilateralFilter->AddObserver( itk::ProgressEvent(), callback );
	bilateralFilter->SetNumberOfThreads( m_threads );
	bilateralFilter->SetRangeSigma(BilatRangeSigma->value()); 
	bilateralFilter->SetDomainSigma(DomainSigma); 
	
	// Connect the pipeline
	bilateralFilter->SetInput( (*m_VectorData)[m_op1->value()].image );
	
	try{
		bilateralFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, bilateralFilter->GetOutput() );
	bilateralFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::MeanFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	InputImageType::SizeType radius;
	
	radius[0] = (unsigned int)mean_radio0->value();
	radius[1] = (unsigned int)mean_radio1->value();
	radius[2] = (unsigned int)mean_radio2->value();
	
	MeanFilterType::Pointer meanFilter = MeanFilterType::New();
	meanFilter->AddObserver( itk::ProgressEvent(), callback );
	meanFilter->SetNumberOfThreads( m_threads );
	meanFilter->SetRadius( radius ); 
	
	// Connect the pipeline
	meanFilter->SetInput( (*m_VectorData)[m_op1->value()].image );
	
	try{
		meanFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, meanFilter->GetOutput() );
	meanFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::WienerFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	WienerFilterType::Pointer wienerFilter = WienerFilterType::New();
	wienerFilter->AddObserver( itk::ProgressEvent(), callback );
	wienerFilter->SetNumberOfThreads( m_threads );
	wienerFilter->Setlambda(lambda->value());
	wienerFilter->Setiterations(static_cast<unsigned int>(iterations->value()));
	wienerFilter->SetNoise(bias->value());
	wienerFilter->SetBias(noise->value());
	wienerFilter->Setgamma(gamma->value());
	wienerFilter->Setopcion(m_filtroWiener-1);
	
	// Connect the pipeline
	wienerFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		wienerFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, wienerFilter->GetOutput() );
	wienerFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::ASRFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	if ( (*m_VectorData)[m_op1->value()].image->GetLargestPossibleRegion().GetSize()[2] == 1){
		ASRFilterType2D::Pointer ASRFilter2D = ASRFilterType2D::New();
		ASRFilter2D->AddObserver( itk::ProgressEvent(), callback );
		ASRFilter2D->SetNumberOfThreads( m_threads );
		ASRFilter2D->SetIter(static_cast<unsigned int>(asr_iterations->value()));
		ASRFilter2D->SetTimeStep( asr_tau->value() );       
		ASRFilter2D->SetBeta( asr_beta->value());
		ASRFilter2D->SetSigma( asr_sigma->value() );
		ASRFilter2D->SetDifference( asr_s->value() );
		
		Cast3Dto2DFilterType::Pointer TDto2DCaster = Cast3Dto2DFilterType::New();
		TDto2DCaster->SetNumberOfThreads( m_threads );
		TDto2DCaster->SetInput( (*m_VectorData)[m_op1->value()].image );
		ASRFilter2D->SetInput( TDto2DCaster->GetOutput() );
		
		Cast2Dto3DFilterType::Pointer DDto3DCaster = Cast2Dto3DFilterType::New();
		DDto3DCaster->SetNumberOfThreads( m_threads );
		DDto3DCaster->SetInput( ASRFilter2D->GetOutput() );
		
		// Update mini-pipeline
		try{
			DDto3DCaster->Update();
		}
		catch( itk::ExceptionObject &e ){
			fl_alert( e.GetDescription() );
			return;
		}
		
		m_VectorData->copyData( m_destino->value()-1, DDto3DCaster->GetOutput() );
		DDto3DCaster = NULL;
	}
	else{
		ASRFilterType::Pointer ASRFilter = ASRFilterType::New();
		ASRFilter->AddObserver( itk::ProgressEvent(), callback );
		ASRFilter->SetNumberOfThreads( m_threads );
		ASRFilter->SetIter(static_cast<unsigned int>(asr_iterations->value()));
		ASRFilter->SetTimeStep( asr_tau->value() );       
		ASRFilter->SetBeta( asr_beta->value());
		ASRFilter->SetSigma( asr_sigma->value() );
		ASRFilter->SetDifference( asr_s->value() );
		
		ASRFilter->SetInput((*m_VectorData)[m_op1->value()].image);
		
		try{
			ASRFilter->Update();
		}
		catch(itk::ExceptionObject &e){
			fl_alert( e.GetDescription() );
			return;
		}
		
		m_VectorData->copyData( m_destino->value()-1, ASRFilter->GetOutput() );
		ASRFilter = NULL;
	}
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}

/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::SRADFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	SRADFilterType::Pointer SRADFilter = SRADFilterType::New();
	SRADFilter->AddObserver( itk::ProgressEvent(), callback );
	SRADFilter->SetNumberOfThreads( m_threads );
	SRADFilter->SetIter(static_cast<unsigned int>(srad_iterations->value()));
	SRADFilter->SetTimeStep( srad_tau->value() );       
	SRADFilter->SetBeta( srad_beta->value());
	SRADFilter->SetExpConstant( srad_expconstant->value() );
	SRADDiffusionType::Pointer SRADdiffusion = SRADDiffusionType::New();
	
	SRADFilter->SetTensorFilter( SRADdiffusion );
	
	// Connect the pipeline
	SRADFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		SRADFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, SRADFilter->GetOutput() );
	SRADFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::DPADFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	DPADFilterType::Pointer DPADFilter = DPADFilterType::New();
	DPADFilter->AddObserver( itk::ProgressEvent(), callback );
	DPADFilter->SetNumberOfThreads( m_threads );
	DPADFilter->SetIter(static_cast<unsigned int>(dpad_iterations->value()));
	DPADFilter->SetTimeStep( dpad_tau->value() );       
	
	// Connect the pipeline
	DPADFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		DPADFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}

	m_VectorData->copyData( m_destino->value()-1, DPADFilter->GetOutput() );
	DPADFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::CurvatureAnisotropicFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	if ( (*m_VectorData)[m_op1->value()].image->GetLargestPossibleRegion().GetSize()[2] == 1){
		Curvature2DDiffusionType::Pointer curvature2D  = Curvature2DDiffusionType::New();
		curvature2D->AddObserver( itk::ProgressEvent(), callback );
		curvature2D->SetNumberOfThreads( m_threads );
		curvature2D->SetNumberOfIterations(static_cast<unsigned int>(AnDifIterations->value()));
		curvature2D->SetConductanceParameter(AnDifConductance->value());
		curvature2D->SetTimeStep(AnDifTimeStep->value());
		
		Cast3Dto2DFilterType::Pointer TDto2DCaster = Cast3Dto2DFilterType::New();
		TDto2DCaster->SetNumberOfThreads( m_threads );
		TDto2DCaster->SetInput( (*m_VectorData)[m_op1->value()].image );
		curvature2D->SetInput( TDto2DCaster->GetOutput() );
		
		Cast2Dto3DFilterType::Pointer DDto3DCaster = Cast2Dto3DFilterType::New();
		DDto3DCaster->SetNumberOfThreads( m_threads );
		DDto3DCaster->SetInput( curvature2D->GetOutput() );
		
		// Update the mini-pipeline:
		try{
			DDto3DCaster->Update();
		}
		catch(itk::ExceptionObject &e){
			fl_alert( e.GetDescription() );
			return;
		}
		
		m_VectorData->copyData( m_destino->value()-1, DDto3DCaster->GetOutput() );
		DDto3DCaster = NULL;
	}
	else{
		CurvatureDiffusionType::Pointer curvature = CurvatureDiffusionType::New();
		curvature->AddObserver( itk::ProgressEvent(), callback );
		curvature->SetNumberOfThreads( m_threads );
		curvature->SetNumberOfIterations( static_cast<unsigned int>(AnDifIterations->value()) );
		curvature->SetConductanceParameter(AnDifConductance->value());
		curvature->SetTimeStep(AnDifTimeStep->value());
		
		curvature->SetInput((*m_VectorData)[m_op1->value()].image);
		
		try{
			curvature->Update();
		}
		catch(itk::ExceptionObject &e){
			fl_alert( e.GetDescription() );
			return;
		}
	
		m_VectorData->copyData( m_destino->value()-1, curvature->GetOutput() );
		curvature = NULL;
	}
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::GradientDiffusionFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	GradientDiffusionType::Pointer gradientDiffusionFilter = GradientDiffusionType::New();
	gradientDiffusionFilter->AddObserver( itk::ProgressEvent(), callback );
	gradientDiffusionFilter->SetNumberOfThreads( m_threads );
	gradientDiffusionFilter->SetNumberOfIterations(static_cast<unsigned int>(GradDifIterations->value()));
	gradientDiffusionFilter->SetConductanceParameter(GradDifConductance->value());
	gradientDiffusionFilter->SetTimeStep(GradDifTimeStep->value());
	gradientDiffusionFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		gradientDiffusionFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, gradientDiffusionFilter->GetOutput() );
	gradientDiffusionFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::GaussianFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();
	gaussianFilter->SetNumberOfThreads( m_threads );
	gaussianFilter->AddObserver( itk::ProgressEvent(), callback );
	gaussianFilter->SetSigma(sigma->value());
	
	// Connect the pipeline
	gaussianFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		gaussianFilter->Update();
	}
	catch( itk::ExceptionObject &e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, gaussianFilter->GetOutput() );
	gaussianFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::MedianFilter( void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	InputImageType::SizeType radius;
	radius[0] = (unsigned int)radio0->value();
	radius[1] = (unsigned int)radio1->value();
	radius[2] = (unsigned int)radio2->value();
	
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();
	medianFilter->AddObserver( itk::ProgressEvent(), callback );
	medianFilter->SetNumberOfThreads( m_threads );
	medianFilter->SetRadius(radius);
	
	// Connect the pipeline
	medianFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		medianFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, medianFilter->GetOutput() );
	medianFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::ABSFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	ABSFilterType::Pointer absFilter = ABSFilterType::New();
	absFilter->AddObserver( itk::ProgressEvent(), callback );
	absFilter->SetNumberOfThreads( m_threads );
	absFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		absFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, absFilter->GetOutput() );
	absFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::ABSValDifFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	ABSValDifFilterType::Pointer absValDifFilter = ABSValDifFilterType::New();
	absValDifFilter->AddObserver( itk::ProgressEvent(), callback );
	absValDifFilter->SetNumberOfThreads( m_threads );
	absValDifFilter->SetInput1((*m_VectorData)[m_op1->value()].image);
	absValDifFilter->SetInput2((*m_VectorData)[m_op2->value()].image);
	
	try{
		absValDifFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, absValDifFilter->GetOutput() );
	absValDifFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::ZeroEdgeFilter(void ){
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	ZeroEdgeFilterType::Pointer zeroEdgeFilter = ZeroEdgeFilterType::New();
	zeroEdgeFilter->AddObserver( itk::ProgressEvent(), callback );
	zeroEdgeFilter->SetNumberOfThreads( m_threads );
	
	// Connect the pipeline
	zeroEdgeFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	//Add window with these two parameters
	//if((z_noparam->value())==0){
	double Zerror = m_max_error->value();
	if ( (Zerror<1) && (Zerror>0) ){
		zeroEdgeFilter->SetMaximumError( m_max_error->value() );
	}
	else{
		if ( Zerror>=1 ){   zeroEdgeFilter->SetMaximumError( 0.9999 );   }
		if ( Zerror<=0 ){   zeroEdgeFilter->SetMaximumError( 0.0001 );   }
    }
	zeroEdgeFilter->SetVariance( m_var->value() );
	//}
	try{
		zeroEdgeFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, zeroEdgeFilter->GetOutput() );
	zeroEdgeFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void FilteringConsole::CannyFilter(void)
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction(this,&FilteringConsole::OnProgress);
	
	CannyFilterType::Pointer cannyFilter = CannyFilterType::New();
	cannyFilter->AddObserver( itk::ProgressEvent(), callback );
	cannyFilter->SetNumberOfThreads( m_threads );
	cannyFilter->SetUpperThreshold( upper_threshold->value() );
	cannyFilter->SetLowerThreshold( lower_threshold->value() );
	
	// Connect the pipeline
	cannyFilter->SetInput( (*m_VectorData)[m_op1->value()].image );
	
	try{
		cannyFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, cannyFilter->GetOutput() );
	cannyFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


	



