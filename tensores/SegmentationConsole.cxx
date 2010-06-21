/*=========================================================================

  Program:   SegmentationConsole.cxx
  Language:  C++
  Date:      27-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#include "SegmentationGUI.h"
#include "SegmentationConsole.h"

// itk includes 
#include <itkImageRegionIteratorWithIndex.h>
//#include <itkGradientImageFilter.h>
#include <itkVTKImageExport.h>
// vtk includes 
#include <vtkITKUtility.h>
#include <vtkImageImport.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
// std includes
#include <vector>

SegmentationConsole::SegmentationConsole(int X, int Y, int W, int H, const char *label)
:SegmentationGUI(X,Y,W,H,label) 
{
  m_thresholder  = BinaryThresholdType::New();
  m_LevelSetSegmentation = SegmentationFilterType::New();
  //m_LevelSetSegmentation->AddObserver( itk::IterationEvent(), iterationCommand);
  m_VTKexporter = VTKexportType::New();
  m_SegmentedCaster = Cast2Dto3DFilterType::New();
  m_InputCaster = Cast3Dto2DFilterType::New();
  m_OutputCaster = Cast2Dto3DFilterType::New();
  m_FloatInputImage = FloatImageType2D::New();
  m_fastMarching = FastMarchingFilter2DType::New();
  
  m_vtkMarcacionElipsoide = vtkMarcacionElipsoide::New();
  m_vtkPlantillaAjustada = vtkPlantillaAjustada::New(); 
  m_vtkFuncionVerosimilitud = vtkFuncionVerosimilitud::New();
  m_vtkOptimizaContorno = vtkOptimizaContorno::New();
  m_vtkValidacion = vtkValidacion::New();
  
  // Para actualizar resultado de Level Set en cada iteracion 
  typedef itk::SimpleMemberCommand< SegmentationConsole > SimpleCommandType;
  SimpleCommandType::Pointer iterationCommand = SimpleCommandType::New();
  iterationCommand->SetCallbackFunction( this, & SegmentationConsole::UpdateViewerAfterIteration );

  // Parametros por defecto de Level Set
  RMSError->value(0.02);
  curvature->value(1);
  propagation->value(1);
  edge->value(1);
  isovalue->value(1);
  lowerThreshold->value(0);
  upperThreshold->value(100);
  reverseExpansionDir->value(0);
  maxCurvatureTimeStep->value(0.25);
  maxPropagationTimeStep->value(0.25);
  diffusionIterations->value(0);
  diffusionConductance->value(1);
  diffusionTimeStep->value(0.125);
  threshIterations->value(0);
  threshConductance->value(0.5);
  threshTimeStep->value(0.1);

}

SegmentationConsole::~SegmentationConsole()
{
}

void SegmentationConsole::SetModelDataBrowser(Fl_Browser *modeldataBrowserIn)
{
	m_modeldataBrowser = modeldataBrowserIn;
}

void SegmentationConsole::SetImageViewer3D(Viewer3D* ImageViewer3DIn)
{
	ImageViewer3D = ImageViewer3DIn;
}

void SegmentationConsole::SetVectorModelData(void* VectorModelData)
{
	m_VectorModelData = (VectorOfModelType*)VectorModelData;
}

void SegmentationConsole::SetVectorData(void* VectorData)
{
	m_VectorData = (VectorOfDataType*)VectorData;
}

void SegmentationConsole::SetImageColorViewer(MyColorViewerType* ImageColorViewerIn)
{
	ImageColorViewer = ImageColorViewerIn;
}

void SegmentationConsole::SetSliders(Fl_Value_Slider **sliders)
{
	sliceNumberSlider = sliders;
}

void SegmentationConsole::SetProgressSliders(fltk::ProgressBar *progressbar)
{
	progressSlider = progressbar;
}

void SegmentationConsole::SetProgressCounter(Fl_Value_Output *progressCounter)
{
	m_ProgressCounter = progressCounter;
}

void SegmentationConsole::SetThreads(int threads)
{
	m_threads = threads;
}

void SegmentationConsole::SetNImagesCargadas(unsigned int &NImagesCargadas)
{
	m_NImagesCargadas = NImagesCargadas;
}

void SegmentationConsole::ClickSelectCallback(float x, float y, float z, float v, int c, void * args )
{
	SegmentationConsole* self = static_cast<SegmentationConsole *>( args );
	self->GetClickPoints(x, y);
}

void SegmentationConsole::GetClickPoints(float x, float y)
{
	seedX->value(x);
	seedY->value(y);
	m_LevelSetSegmentation->Modified();
}

void SegmentationConsole::OnProgress(itk::Object *object, const itk::EventObject &event)
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

/* Modified */
void SegmentationConsole::LoadME( void ){
	std::cout << "En LoadME" << std::endl;
	ME_filename = fl_file_chooser("Image Filename","*.{vtk,txt}","");
	ME_FilePageBrowser->value(ME_filename);
}

/* Modified */
void SegmentationConsole::LoadKNN( void ){
	std::cout << "En LoadKNN" << std::endl;
	KNN_filename = fl_file_chooser("Points Filename","*.{txt,mhd}","");
	KNN_FilePageBrowser->value(KNN_filename);
}

/* Modified */
void SegmentationConsole::KNNFileInputChange() {
	const char *text = KNN_FilePageBrowser->value();
	if (text != NULL && strlen(text) > 0){
		KNN_filename = text;
	}
}


/* Modified */
void SegmentationConsole::KnnFilter(void)
{
	float x,y,z;
	std::cout << "Haciendo kNN " << std::endl;
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// Make use of the progress bar!!!
	//  CommandPointer callback = CommandType::New();
	//  callback->SetCallbackFunction(this,&UsimagToolConsole::OnProgress);
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
	IndexType v_index;
	v_index.Fill(0);
	
	VoronoiInputType::SizeType   size;
	VoronoiInputType::RegionType diagramRegion;
	diagramRegion.SetIndex(v_index);
	
	//Minimum of the input image:
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// Move this definition to an appropriate location!!!
	typedef itk::MinimumMaximumImageCalculator<InputImageType> MinMaxCalculatorType;
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	MinMaxCalculatorType::Pointer Ch1Calculator = MinMaxCalculatorType::New();
	Ch1Calculator->SetImage(  (*m_VectorData)[m_op1->value()].image  );
	try{
		Ch1Calculator->ComputeMaximum();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	size[0] = (unsigned int)Ch1Calculator->GetMaximum()+1;
	if( size[0]<256 ){size[0]=256;}
	diagramRegion.SetSize(size);
	std::cout << "inputRegion del Test = " << size << std::endl;
	
	VoronoiInputType::Pointer inputvoronoi = VoronoiInputType::New();
	inputvoronoi->SetRegions(diagramRegion);
	inputvoronoi->Allocate();
	inputvoronoi->FillBuffer(-1);
	
	classifyKNNCoreType::ListType proto;
	classifyKNNCoreType::NodeType auxNode;
	DiagramImageType::IndexType   index;
	
	char buffer[2048];
	int pclass;
	
	FILE* fp = fopen( (const char*)KNN_filename, "r" );
	while( !feof(fp) ){
		fgets( buffer, sizeof(buffer), fp );
		float ind_aux;
		sscanf( buffer, "%f %f %f %f %d", &x, &y, &z, &ind_aux, &pclass );
		index[0]                = (int)ind_aux;
		auxNode.index_intensity = index;
		auxNode.clase           = pclass;
		proto.push_back(auxNode);
	}
	fclose(fp);
	
	//std::cout << "proto.size() " << proto.size() << std::endl;
	for( unsigned int i=0; i<proto.size(); i++ )
		inputvoronoi->SetPixel( proto[i].index_intensity, i );
	
	//std::cout << "Proto[0], clase:" << proto[0].clase << " index " << proto[0].index_intensity << std::endl;
	VoronoiFilterType::Pointer voronoi = VoronoiFilterType::New();
	voronoi->SetBackgroundValue( -1 );
	voronoi->SetInput( inputvoronoi );
	voronoi->SetK(   static_cast<unsigned int>( KNN_K->value() )   );
	try{
		voronoi->Update();
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	classifyKNNCoreType::Pointer knnFilter = classifyKNNCoreType::New();
	knnFilter->SetNumberOfThreads( m_threads );
	knnFilter->SetDiagram( voronoi->GetOutput() );
	knnFilter->SetCh1(  (*m_VectorData)[m_op1->value()].image  );
	knnFilter->SetPrototipos( proto );
	knnFilter->SetK(   static_cast<unsigned int>( KNN_K->value() )   );
	try{
		knnFilter->Update();;
	}
	catch(itk::ExceptionObject &e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, knnFilter->GetOutput() );
	knnFilter = NULL;
	
	// Reset progress bar
	progressSlider->value(0);
	m_ProgressCounter->value(0);
	std::cout << "KNN Done" << std::endl;
}

void SegmentationConsole::MarcacionElipseInit()
{
	std::cout << "Doing MarcacionElipseInit " << ME_filename << std::endl;
	m_vtkPlantillaAjustada->SetFichModelo(ME_filename);
	m_vtkFuncionVerosimilitud->SetRuta(ME_filename);
	m_vtkMarcacionElipsoide->SetNomFichMar(ME_filename);
	std::cout << "Setting values" << std::endl;
	m_vtkPlantillaAjustada->SetJ( static_cast<int>(ME_J->value()));
	m_vtkFuncionVerosimilitud->SetJ( static_cast<int>(ME_J->value()));
	m_vtkPlantillaAjustada->SetK( static_cast<int>(ME_K->value()));
	m_vtkFuncionVerosimilitud->SetK( static_cast<int>(ME_K->value()));
	m_vtkPlantillaAjustada->Setdrmax( static_cast<float>(ME_drmax->value()));
	m_vtkFuncionVerosimilitud->Setdrmax( static_cast<float>(ME_drmax->value()));
	m_vtkFuncionVerosimilitud->SetNg(static_cast<int>(ME_Ng->value()));
	m_vtkPlantillaAjustada->SetUsimagTool(1);
	m_vtkFuncionVerosimilitud->SetUsimagTool(1);
	m_vtkOptimizaContorno->SetUsimagTool(1);
}

void SegmentationConsole::CalculoModeloVerosimilitud() {

  double m_spac_in[3],m_spac_out[3];
  int i,j;
   
  std::cout << "Doing CalculoModeloVerosimilitud" << std::endl;
   //Hay que incluir vtkPolyData.h, vtkMatrix4x4.h, vtkPoints.h, vtkImageData.h, vtkFloatArray.h (... y quiza alguno mas que se me haya pasado).
   
  vtkPolyData* poli = vtkPolyData::New();

  if (static_cast<int>(ME_Carga->value()) == 1) {
    poli = m_vtkMarcacionElipsoide->GeneraElipsoide( m_vtkMarcacionElipsoide->CargarMarcacionInicial());
  } else {
    // poli = m_vtkMarcacionElipsoide->GeneraElipsoide( ME_Fiducials );
    //m_vtkMarcacionElipsoide->GuardarMarcacionInicial( ME_Fiducials );
  }
   
   vtkImageImport* vtkImporter = vtkImageImport::New();
  m_VTKexporter->SetInput((*m_VectorData)[m_op1->value()].image);
  ConnectPipelines(m_VTKexporter, vtkImporter);
  vtkImporter->GetOutput()->Update();

  vtkImporter->GetOutput()->GetSpacing(m_spac_in);
   
  for (i=0;i<3;i++)
     m_spac_out[(i+m_activeinput)%3]=m_spac_in[i]*1000; //Cambio de m. a mm. y correspondencia dimensional.
  
  
  vtkMatrix4x4 *plan=vtkMatrix4x4::New();
  vtkMatrix4x4 *tran=vtkMatrix4x4::New();
  for (i=0;i<4;i++)  {
    for (j=0;j<4;j++) {
      plan->SetElement(i,j,0);
      tran->SetElement(i,j,0);
    }
  }
  plan->SetElement(1,0,-1);
  plan->SetElement(2,1,1);
  plan->SetElement(0,2,-1);
  plan->SetElement(3,3,1);
  tran->SetElement(0,1,-1/m_spac_out[0]);
   tran->SetElement(1,2,1/m_spac_out[1]);
  tran->SetElement(2,0,1/m_spac_out[2]);
  for (i=0;i<3;i++)
     tran->SetElement(i,3,(int)(*m_VectorData)[m_op1->value()].image->GetLargestPossibleRegion().GetSize()[i]/2);
  tran->SetElement(3,3,1);
   
  //std::cout << "slide number " << ImageViewer[m_activeinput]->GetNumSlices() << std::endl;

  m_vtkPlantillaAjustada->IntroducePlano(plan);
  m_vtkPlantillaAjustada->SetSlth(m_spac_out[2]);
  m_vtkPlantillaAjustada->Setdrmax(m_vtkPlantillaAjustada->GetSlth()*m_vtkPlantillaAjustada->Getdrmax());
  m_vtkPlantillaAjustada->ObtieneParametros(poli);
  

  m_vtkFuncionVerosimilitud->SetP(m_vtkPlantillaAjustada->GetP());
  m_vtkFuncionVerosimilitud->Setdrmax(m_vtkPlantillaAjustada->Getdrmax());
  m_vtkFuncionVerosimilitud->IntroduceMatriza(tran);
  m_vtkFuncionVerosimilitud->IntroducePlano(plan);
  m_vtkFuncionVerosimilitud->Setlimite(m_vtkPlantillaAjustada->Getlimite());
   
  m_vtkFuncionVerosimilitud->IntroduceImagen( vtkImporter->GetOutput(), 0);
  m_vtkFuncionVerosimilitud->IntroduceParam(m_vtkPlantillaAjustada->ObtieneParam());
  m_vtkFuncionVerosimilitud->IntroduceCentro(m_vtkPlantillaAjustada->ObtieneCentro());
  m_vtkFuncionVerosimilitud->Ejecuta();
	
  m_vtkPlantillaAjustada->EstableceMascaraElipsoide(m_vtkFuncionVerosimilitud->DevuelveMascaraElipsoide());
  
  // Renderizar!!!!

}

void SegmentationConsole::CalculoOptimizacionSA() {
   
   //Quiza incluir vtkStructuredPoints, vtkPolyData, vtkPolyDataWriter
   
   //Ahi va una posibilidad para guardar en fichero:
   vtkPolyData *Poli=vtkPolyData::New();
   vtkPolyDataWriter *PolWri=vtkPolyDataWriter::New();

    std::cout << "Doing OptimizacionSA" << std::endl;
   
  m_vtkOptimizaContorno->SetSuper(0);
  m_vtkOptimizaContorno->SetPeriodMalla(1,0,0);
  m_vtkOptimizaContorno->SetNumEntidades(1);
  m_vtkOptimizaContorno->SetDimensionalidadEstados(1);
  m_vtkOptimizaContorno->SetDimensionalidadMalla(2);
  m_vtkOptimizaContorno->SetIndependencia(0);
  m_vtkOptimizaContorno->SetOrdenMalla(4);
  m_vtkOptimizaContorno->SetOrdenEstados(5);
  m_vtkOptimizaContorno->SetDimMalla(m_vtkFuncionVerosimilitud->GetJ(),m_vtkFuncionVerosimilitud->GetP(),1);

  m_vtkOptimizaContorno->SetDimZ(1,1,1);
  m_vtkOptimizaContorno->SetK(m_vtkFuncionVerosimilitud->GetK());
  m_vtkOptimizaContorno->SetV(2); 
  m_vtkOptimizaContorno->ConstruyeVecindario();
  m_vtkOptimizaContorno->EstableceLE(m_vtkFuncionVerosimilitud->DevuelveLR(),0,0);
	

  m_vtkOptimizaContorno->Llamada(0);
 std::cout << "Finished OptimizacionSA" << std::endl;   
   //Prosigue el almacenamiento en fichero.
   
   Poli=m_vtkPlantillaAjustada->ConstruyeModelo(m_vtkOptimizaContorno->DevuelveRho());
   Poli->Update();
   
   //Todo esto en un if con variable para guardar en fichero.
   PolWri->SetFileName("prueba_resultados.vtk");
   PolWri->SetFileTypeToASCII();
   PolWri->SetInput(Poli);
   PolWri->Write();
   
    std::cout << "Finished FileWriting" << std::endl;
   
   //Renderizacion:
   
  //set MarcacionElipsoide(Poli,$MarcacionElipsoide(color)) [MarcacionElipsoide(Plantilla,$id) ConstruyeModelo [MarcacionElipsoide(Optimiza,$id) DevuelveRho]]
  //set max [MarcacionElipsoide(Plantilla,$id) GetMax]
  //set min [MarcacionElipsoide(Plantilla,$id) GetMin]
  //MarcacionElipsoide(Generador,$id) SetMax [lindex $max 0] [lindex $max 1] [lindex $max 2] [lindex $max 3]
  //MarcacionElipsoide(Generador,$id) SetMin [lindex $min 0] [lindex $min 1] [lindex $min 2] [lindex $min 3]
  //set tipo 0
  //set MarcacionElipsoide(actor,$MarcacionElipsoide(color)) [MarcacionElipsoide(Generador,$id) DibujaModelo $MarcacionElipsoide(Poli,$MarcacionElipsoide(color)) $tipo $MarcacionElipsoide(color)]	
  //MainAddActor $MarcacionElipsoide(actor,$MarcacionElipsoide(color))

  //if { $tipo==1 } { 
  //set MarcacionElipsoide(actor,[expr 1+$MarcacionElipsoide(color)]) [MarcacionElipsoide(Generador,$id) DevuelveBarra]
  //MainAddActor $MarcacionElipsoide(actor,[expr 1+$MarcacionElipsoide(color)]) }

  //RenderAll
  
}

/* Modified */
void SegmentationConsole::ThresholdFilter( void )
{
	// Create the segmentation filter
	CommandPointer callback = CommandType::New();  
	callback->SetCallbackFunction(this,&SegmentationConsole::OnProgress);
	
	BinaryThresholdFilterType::Pointer thresholdFilter = BinaryThresholdFilterType::New();
	thresholdFilter->AddObserver( itk::ProgressEvent(), callback );
	thresholdFilter->SetNumberOfThreads( m_threads );
	thresholdFilter->SetUpperThreshold( threshold_above->value() ); 
	thresholdFilter->SetLowerThreshold( threshold_below->value() ); 
	thresholdFilter->SetInsideValue( 1 );
	thresholdFilter->SetOutsideValue( 0 );
	thresholdFilter->SetInput(   (*m_VectorData)[m_op1->value()].image   );
	
	try{
		thresholdFilter->Update();
	}
	catch( itk::ExceptionObject &e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, thresholdFilter->GetOutput() );
	thresholdFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}



/* Modified */
void SegmentationConsole::OtsuThresholdFilter( void )
{
	// Create the segmentation filter
	CommandPointer callback = CommandType::New();  
	callback->SetCallbackFunction(this,&SegmentationConsole::OnProgress);
	
	OtsuThresholdFilterType::Pointer OtsuthresholdFilter = OtsuThresholdFilterType::New();
	OtsuthresholdFilter->AddObserver( itk::ProgressEvent(), callback );
	OtsuthresholdFilter->SetNumberOfThreads( m_threads );
	OtsuthresholdFilter->SetNumberOfHistogramBins( bins->value() ); 
	OtsuthresholdFilter->SetInsideValue( 1 );
	OtsuthresholdFilter->SetOutsideValue( 0 );
	OtsuthresholdFilter->SetInput(   (*m_VectorData)[m_op1->value()].image   );
	
	try{
		OtsuthresholdFilter->Update();
	}
	catch( itk::ExceptionObject &e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_VectorData->copyData( m_destino->value()-1, OtsuthresholdFilter->GetOutput() );
	OtsuthresholdFilter = NULL;
	
	progressSlider->value(0);
	m_ProgressCounter->value(0);
}


void SegmentationConsole::ResetAllParameters() {

  m_thresholder->SetOutsideValue(0);
  m_thresholder->SetInsideValue(1);
  m_thresholder->SetUpperThreshold(0.5);
  m_thresholder->SetLowerThreshold(-0.5);
  
  //m_Curvature->SetNumberOfIterations(static_cast<unsigned int>(diffusionIterations->value()));
  //m_Curvature->SetConductanceParameter(diffusionConductance->value());
  //m_Curvature->SetTimeStep(diffusionTimeStep->value());
  
  m_LevelSetSegmentation->AbortGenerateDataOff();
  m_LevelSetSegmentation->SetUpperThreshold( upperThreshold->value() );
  m_LevelSetSegmentation->SetLowerThreshold( lowerThreshold->value() );
  m_LevelSetSegmentation->SetCurvatureScaling( curvature->value() );
  m_LevelSetSegmentation->SetPropagationScaling( propagation->value() );
  m_LevelSetSegmentation->SetEdgeWeight( edge->value() );
  m_LevelSetSegmentation->SetMaximumRMSError( RMSError->value() );
  m_LevelSetSegmentation->SetNumberOfIterations( (unsigned int)levelsetiterations->value() );
  m_LevelSetSegmentation->SetMaximumCurvatureTimeStep( maxCurvatureTimeStep->value() );
  m_LevelSetSegmentation->SetMaximumPropagationTimeStep( maxPropagationTimeStep->value() );
  m_LevelSetSegmentation->SetSmoothingIterations( (unsigned int)threshIterations->value() );
  m_LevelSetSegmentation->SetSmoothingConductance( threshConductance->value() );
  m_LevelSetSegmentation->SetSmoothingTimeStep( threshTimeStep->value() );

  if(reverseExpansionDir->value()) {
    m_LevelSetSegmentation->ReverseExpansionDirectionOn();
  }
  else {
    m_LevelSetSegmentation->ReverseExpansionDirectionOff();
  }
}

void SegmentationConsole::UpdateViewerAfterIteration(){
  static unsigned int iterationCounter = 0;
  //std::cout << "UpdateViewerAfterIteration  " << "it n:" <<  m_LevelSetSegmentation->GetElapsedIterations() << std::endl;
  if( (iterationCounter%((int)updateIterations->value()) == 0)
      && (iterationCounter != 0) ) {
    // Move the pixel container and image information of the image we are working
    // on into a temporary image to  use as the input to the mini-pipeline.  This
    // avoids a complete copy of the image.
    FloatImageType2D::Pointer tmp = FloatImageType2D::New();
    tmp->SetRequestedRegion( m_LevelSetSegmentation->GetOutput()->GetRequestedRegion() );
    tmp->SetBufferedRegion( m_LevelSetSegmentation->GetOutput()->GetBufferedRegion() );
    tmp->SetLargestPossibleRegion( m_LevelSetSegmentation->GetOutput()->GetLargestPossibleRegion() );
    tmp->SetPixelContainer( m_LevelSetSegmentation->GetOutput()->GetPixelContainer() );
    tmp->CopyInformation( m_LevelSetSegmentation->GetOutput() );

    // update overlay
    m_thresholder->Modified();
    m_thresholder->SetInput(tmp);
    try
      {
      // Updates m_thresholder
      m_SegmentedCaster->UpdateLargestPossibleRegion();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      }

    // update the viewer with the current segmentation
    //ImageViewer[m_activeinput]->SetOverlay(m_SegmentedCaster->GetOutput());
    //ImageViewer[m_activeinput]->SetOverlayOpacity(1.0);
    //ImageViewer[m_activeinput]->Update();
    ImageViewer[m_activeinput]->SetOverlay(m_SegmentedCaster->GetOutput());
    ImageViewer[m_activeinput]->SetOverlayOpacity(1.0);
    ImageViewer[m_activeinput]->Update();
    Fl::check();

    elapsedIterations->value(m_LevelSetSegmentation->GetElapsedIterations());
    lastRMSChange->value(m_LevelSetSegmentation->GetRMSChange());
    progressSlider->value( m_LevelSetSegmentation->GetElapsedIterations()
			   / levelsetiterations->value() );
    m_ProgressCounter->value(100*m_LevelSetSegmentation->GetElapsedIterations()
			     / levelsetiterations->value());
  }
  iterationCounter++;
}

void SegmentationConsole::ThresholdLevelSetSegmentation(void ) {
  std::cout << "Doing Level Set 2D" << std::endl;
 
  // CommandPointer callback = CommandType::New();  
  // callback->SetCallbackFunction(this,&UsimagToolConsole::OnProgress);
  //m_LevelSetSegmentation->AddObserver( itk::ProgressEvent(), callback);

  this->ResetAllParameters();

  m_InputCaster->SetInput((*m_VectorData)[m_op1->value()].image);
  m_InputCaster->Update();
  m_FloatInputImage = m_InputCaster->GetOutput();
  m_FloatInputImage->SetSpacing(1.0);

  m_LevelSetSegmentation->SetFeatureImage( m_FloatInputImage);
  m_LevelSetSegmentation->SetInput( m_fastMarching->GetOutput() );
  m_LevelSetSegmentation->SetIsoSurfaceValue(0.0);

  m_seedPosition2D[0] = (unsigned long)seedX->value();
  m_seedPosition2D[1] = (unsigned long)seedY->value();
  
  m_node2D.SetValue( -distance->value() );
  m_node2D.SetIndex( m_seedPosition2D );
  
  m_seeds2D->Initialize();
  m_seeds2D->InsertElement( 0, m_node2D );
  m_fastMarching->SetTrialPoints( m_seeds2D );
 
  m_SegmentedCaster->SetInput(m_thresholder->GetOutput());
  m_fastMarching->SetOutputSize(m_FloatInputImage->GetBufferedRegion().GetSize() );
  
  //std::cout << "Iterations " << levelsetiterations->value() << std::endl;

  //m_FloatImage->SetRegions(m_FloatInputImage->GetLargestPossibleRegion());
  //m_FloatImage->Allocate();
  //m_FloatImage = m_LevelSetSegmentation->GetOutput();
  //m_FloatImage = const_cast<FloatImageType *>(m_LevelSetSegmentation->GetSpeedImage());
  
  try {
    m_LevelSetSegmentation->UpdateLargestPossibleRegion();    
  }
  catch(itk::ExceptionObject &e) {
    std::cerr << "Caught ITK exception: " << e << std::endl;
    std::cerr << e << std::endl;
  }

  //m_real_to_char->SetInput(m_FloatImage);
  //m_real_to_char->Update();
  //m_OutputImage = m_real_to_char->GetOutput();
  
  m_OutputCaster->SetInput(m_LevelSetSegmentation->GetOutput());
  m_OutputCaster->Update();
  
  m_VectorData->copyData( m_destino->value()-1, m_OutputCaster->GetOutput() );

  progressSlider->value(0);
  m_ProgressCounter->value(0);
}
