/*=========================================================================

  Program:   BasicOpConsole.cxx
  Language:  C++
  Date:      7-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#include "BasicOpGUI.h"
#include "BasicOpConsole.h"
// itk includes 
#include <itkImageRegionIteratorWithIndex.h>
#include <itkVTKImageExport.h>
// vtk includes 
#include <vtkITKUtility.h>
#include <vtkImageImport.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
// std includes
#include <vector>

/** Last modified: 13/07/2008 AN-TONIO */
BasicOpConsole::BasicOpConsole(int X, int Y, int W, int H, const char *label)
:BasicOpGUI(X,Y,W,H,label) 
{
}

BasicOpConsole::~BasicOpConsole()
{
}

void BasicOpConsole::SetModelDataBrowser(Fl_Browser *modeldataBrowserIn)
{
	m_modeldataBrowser = modeldataBrowserIn;
}

void BasicOpConsole::SetImageViewer3D(Viewer3D* ImageViewer3DIn)
{
	ImageViewer3D = ImageViewer3DIn;
}

void BasicOpConsole::SetVectorModelData(void* VectorModelData){
	m_VectorModelData = (VectorOfModelType*)VectorModelData;
}

void BasicOpConsole::SetVectorData(void* VectorData){
	m_VectorData = (VectorOfDataType*)VectorData;
}

void BasicOpConsole::SetImageColorViewer(MyColorViewerType* ImageColorViewerIn)
{
	ImageColorViewer = ImageColorViewerIn;
}

void BasicOpConsole::SetThreads(int threads)
{
	m_threads = threads;
}

void BasicOpConsole::SetProgressSliders(fltk::ProgressBar *progressbar)
{
	progressSlider = progressbar;
}

void BasicOpConsole::SetProgressCounter(Fl_Value_Output *progressCounter)
{
	m_ProgressCounter = progressCounter;
}

void BasicOpConsole::SetSliders(Fl_Value_Slider **sliders)
{
	sliceNumberSlider = sliders;
}

void BasicOpConsole::OnProgress(itk::Object *object, const itk::EventObject &event )
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
void BasicOpConsole::AddFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image || !(*m_VectorData)[m_op2->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	AddFilterType::Pointer addFilter = AddFilterType::New();
	addFilter->AddObserver( itk::ProgressEvent(), callback );
	addFilter->SetNumberOfThreads( m_threads );
	addFilter->SetInput1((*m_VectorData)[m_op1->value()].image);
	addFilter->SetInput2((*m_VectorData)[m_op2->value()].image);
	
	try{
		addFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	m_VectorData->copyData( m_destino->value()-1, addFilter->GetOutput() );
	addFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void BasicOpConsole::MultiplyFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image || !(*m_VectorData)[m_op2->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
	multiplyFilter->AddObserver( itk::ProgressEvent(), callback );
	multiplyFilter->SetNumberOfThreads( m_threads );
	multiplyFilter->SetInput1((*m_VectorData)[m_op1->value()].image);
	multiplyFilter->SetInput2((*m_VectorData)[m_op2->value()].image);
	
	try{
		multiplyFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	m_VectorData->copyData( m_destino->value()-1, multiplyFilter->GetOutput() );
	multiplyFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


/** Last modified: 13/07/2008 AN-TONIO */
void BasicOpConsole::RelabelFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	relabelFilter->AddObserver( itk::ProgressEvent(), callback );
	relabelFilter->SetNumberOfThreads( m_threads );
	relabelFilter->SetChange( relabel_orig->value(), relabel_result->value() ); 
	relabelFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		relabelFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, relabelFilter->GetOutput() );
	relabelFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/** Last modified: 13/07/2008 AN-TONIO */
void BasicOpConsole::RescaleFilter(void )
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->AddObserver( itk::ProgressEvent(), callback );  
	rescaleFilter->SetNumberOfThreads( m_threads );
	rescaleFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	rescaleFilter->SetOutputMaximum(rescale_max->value());
	rescaleFilter->SetOutputMinimum(rescale_min->value());
	
	try{
		rescaleFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, rescaleFilter->GetOutput() );
	rescaleFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/** Last modified: 13/07/2008 AN-TONIO */
void BasicOpConsole::ErodeFilter(void ) 
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	StructuringElementType structuringElement;
	InputImageType::SizeType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;
	structuringElement.SetRadius( radius );  // 3x3 ball structuring element
	structuringElement.CreateStructuringElement();
	std::cout << "size: " << structuringElement.Size() << " radius: " << structuringElement.GetRadius() << std::endl;
	
	
	ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
	erodeFilter->AddObserver( itk::ProgressEvent(), callback );
	erodeFilter->SetNumberOfThreads( m_threads );
	// Connect the pipeline
	erodeFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	erodeFilter->SetKernel( structuringElement );
	erodeFilter->SetObjectValue( m_objectvalue->value() );
	
	try{
		erodeFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, erodeFilter->GetOutput() );
	erodeFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/** Last modified: 13/07/2008 AN-TONIO */
void BasicOpConsole::DilateFilter(void ) 
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	StructuringElementType structuringElement;
	structuringElement.SetRadius( 1 );  // 3x3 structuring element
	structuringElement.CreateStructuringElement();
	
	DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
	dilateFilter->AddObserver( itk::ProgressEvent(), callback );  
	dilateFilter->SetNumberOfThreads( m_threads );
	// Connect the pipeline
	dilateFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	dilateFilter->SetKernel( structuringElement );
	dilateFilter->SetObjectValue( m_objectvalue2->value() );
	
	try{
		dilateFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, dilateFilter->GetOutput() );
	dilateFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/** Last modified: 13/07/2008 AN-TONIO */
void BasicOpConsole::RelabelCompFilter(void ) 
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	RelabelCompFilterType::Pointer relabelCompFilter = RelabelCompFilterType::New();
	relabelCompFilter->AddObserver( itk::ProgressEvent(), callback );
	relabelCompFilter->SetNumberOfThreads( m_threads );
	// Connect the pipeline
	relabelCompFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		relabelCompFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, relabelCompFilter->GetOutput() );
	relabelCompFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/** Last modified: 13/07/2008 AN-TONIO */
void BasicOpConsole::InvertFilter(void ) 
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	InvertFilterType::Pointer invertFilter = InvertFilterType::New();
	invertFilter->AddObserver( itk::ProgressEvent(), callback );
	invertFilter->SetNumberOfThreads( m_threads );
	invertFilter->SetNumberOfThreads(m_threads);
	invertFilter->SetMaximum(max_out->value());
	
	// Connect the pipeline
	invertFilter->SetInput((*m_VectorData)[m_op1->value()].image);
	
	try{
		invertFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, invertFilter->GetOutput() );
	invertFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/** Last modified: 13/07/2008 AN-TONIO */
void BasicOpConsole::GradientFilter(void)
{
	if( !(*m_VectorData)[m_op1->value()].image ){return;}
	CommandPointer callback = CommandType::New();
	callback->SetCallbackFunction( this, &BasicOpConsole::OnProgress );
	
	GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
	gradientFilter->SetInput( (*m_VectorData)[m_op1->value()].image );
	gradientFilter->AddObserver( itk::ProgressEvent(), callback );
	gradientFilter->SetNumberOfThreads( m_threads );
	
	try{
		gradientFilter->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, gradientFilter->GetOutput() );
	gradientFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



