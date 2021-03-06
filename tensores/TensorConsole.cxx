/*=========================================================================
 
 Program:   TensorConsole.cxx
 Language:  C++
 Date:      7-05-2008
 Version:   1.0
 
 Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
 See http://www.lpi.tel.uva.es/UsimagTool for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE. 
 
 =========================================================================*/
#include "TensorGUI.h"
#include "TensorConsole.h"
#include "geodesicPath3D.h"
#include "tensor/itkDWImages.h"
// itk includes 
#include <itkImageRegionIteratorWithIndex.h>
#include <itkGradientImageFilter.h>
#include <itkVTKImageExport.h>
#include <itkImageFileWriter.h>
// vtk includes 
#include <vtkITKUtility.h>
#include <vtkImageImport.h>
#include <vtkIdList.h>
#include <vtkRungeKutta4.h>
#include <vtkStreamLine.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkCardinalSpline.h>
#include <vtkPolyLine.h>
#include <vtkActor.h>
#include <vtkProperty.h>
// std includes
#include <vector>

//Añadido por mí. OJO!! No hay que quitarlos TODOS!!
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLight.h>
#include <vtkLightCollection.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkTensorGlyph.h>
#include <vtkStructuredPoints.h>
#include <vtkSphereSource.h>
#include <vtkMath.h>
#include <vtkCubeSource.h>
#include <vtkPolyDataSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkSuperquadricSource.h>
#include <vtkLODActor.h>

//Añadidos por supercuádricas
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkScalarBarActor.h"
#include "vtkLookupTable.h"
#include "vtkProp.h"
#include "vtkColorTransferFunction.h"
#include <math.h>
#include "vtkTensorGlyphDTI.h"

#include "vtkPolyDataMapper2D.h"
#include "vtkActor2D.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkProperty2D.h"

#include "vtkCylinderSource.h"
#include "vtkConeSource.h"
#include "strain/itkStrainTensor.h"

#include "strain/vtkTensorGlyphStrain.h"
#include "vtkGlyphSource2D.h"

#include "time.h"

TensorConsole::TensorConsole(int X, int Y, int W, int H, const char *label)
:TensorGUI(X,Y,W,H,label) 
{
	m_seeds = NodeContainer::New();
	m_geodesicPath3D = geodesicPath3DType::New();
	m_dwifilter = DWIFilterType::New();
	m_dtiEstimator = DTIEstimatorType::New();
	m_HasBeenRegistered=false;
	m_HasBeenFiberTracked=false;
	
	m_register=NULL;
	m_SeedRegionsImage=NULL;
	m_TractRegionsImage=NULL;
	
	// inicializamos los text buffers para medidas en fibras
	for (int i=0;i<19;i++) m_buff[i] = new Fl_Text_Buffer();
	m_buffNfibers = new Fl_Text_Buffer();
	m_bufferROIsText = new Fl_Text_Buffer();

	m_VTKexporter = VTKexportType::New();
	
	rangeStrainMin = 0.0;
	rangeStrainMax = 3.0;

	//==========================================================================
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** This piece of code is to disable the menus and functionalities that physicians do not require */
#ifdef FORPHYSICIANS
	// Disable the menus corresponding to Basic Ops., Filtering, Segmentation and Registration:
	GeoPathButton->deactivate();
#endif
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
	/** *********************************************************************************************************************** */
}

TensorConsole::~TensorConsole() {
	
}

void TensorConsole::SetModelDataBrowser(Fl_Browser *modeldataBrowserIn)
{
	m_modeldataBrowser = modeldataBrowserIn;
}

void TensorConsole::SetDWIDataBrowser(Fl_Browser *DWIdataBrowserIn)
{
	m_DWIdataBrowser = DWIdataBrowserIn;
}

void TensorConsole::SetStrainTensorDataBrowser(Fl_Browser *STDataBrowserIn)
{
	m_strainTensorDataBrowser = STDataBrowserIn;
}

//void TensorConsole::SetImageViewer3D(Viewer3D* ImageViewer3DIn)
//{
//	ImageViewer3D = ImageViewer3DIn;
//}

void TensorConsole::SetVectorModelData(void* VectorModelData)
{
	m_VectorModelData = (VectorOfModelType*)VectorModelData;
}

void TensorConsole::SetVectorTensorData(void* VectorTensorData)
{
	m_VectorTensorData = (VectorOfTensorDataType*)VectorTensorData;
}

void TensorConsole::SetVectorData(void* VectorData)
{
	m_VectorData = (VectorOfDataType*)VectorData;
}

void TensorConsole::SetVectorDWIData(void* VectorDWIData)
{
	m_VectorDWIData = (VectorOfDWIDataType*)VectorDWIData;
}

void TensorConsole::SetVectorSTData(void* VectorSTData)
{
	m_VectorSTData = (VectorOfSTDataType*)VectorSTData;
}

void TensorConsole::SetImageColorViewer(MyColorViewerType* ImageColorViewerIn)
{
	ImageColorViewer = ImageColorViewerIn;
}

void TensorConsole::SetImageOverlay(UCharImageType* m_imageOverlayIn)
{
	m_imageOverlay = m_imageOverlayIn;
}

void TensorConsole::SetSliders(Fl_Value_Slider **sliders){
	sliceNumberSlider = sliders;
}

void TensorConsole::SetPanelMessage(Fl_Double_Window *panelmessageIn) {
	panelMessage = panelmessageIn;
}

void TensorConsole::SetOrientationText(Fl_Text_Display *orientationTextIn) {
	m_OrientationText = orientationTextIn;
}
void TensorConsole::SetOriginText(Fl_Text_Display *originTextIn) {
	m_OriginText = originTextIn;
}
void TensorConsole::SetMinScalarRange(Fl_Value_Slider *min_scalar_range_In) {
	min_scalar_range = min_scalar_range_In;
}
void TensorConsole::SetMaxScalarRange(Fl_Value_Slider *max_scalar_range_In) {
	max_scalar_range = max_scalar_range_In;
}
void TensorConsole::SetActiveStrainImage(vtkImageData *image) {
	m_activeStrainImage = image;
}

void TensorConsole::InfoImagen( ) {  
	int dataId = m_tensordataBrowser->value()-1;
	if (dataId < 0) return;
	char texto[100];
	float v[3][3];
	for (unsigned int i=0;i<3;i++) {
		for (unsigned int j=0;j<3;j++) {
			v[i][j] = (float)(*m_VectorTensorData)[dataId].image->GetDirection()[i][j];
		}
	}
	
	Fl_Text_Buffer *buff = new Fl_Text_Buffer();
	sprintf(texto,"Orientation \n %f, %f, %f \n %f %f %f \n %f %f %f",v[0][0],v[0][1],v[0][2],v[1][0],v[1][1],v[1][2],v[2][0],v[2][1],v[2][2]);
	buff->text(texto);
	m_OrientationText->buffer(buff);
	
	char texto2[100];
	float origen[3];
	origen[0] = (float)(*m_VectorTensorData)[dataId].image->GetOrigin()[0];
	origen[1] = (float)(*m_VectorTensorData)[dataId].image->GetOrigin()[1];
	origen[2] = (float)(*m_VectorTensorData)[dataId].image->GetOrigin()[2];
	Fl_Text_Buffer *buff2 = new Fl_Text_Buffer();
	sprintf(texto2,"Origin %f, %f, %f",origen[0],origen[1],origen[2]);
	buff2->text(texto2);
	m_OriginText->buffer(buff2);
	
	//char texto3[100];
	//	typedef itk::MinimumMaximumImageCalculator< InputImageType > MinMaxCalculatorType;
	//	MinMaxCalculatorType::Pointer calculator = MinMaxCalculatorType::New( );
	//	calculator->SetImage( m_VectorTensorData[dataId].image );
	//	
	//	calculator->Compute( );
	//	float minIntensity = calculator->GetMinimum( );
	//	float maxIntensity = calculator->GetMaximum( );
	//	std::cout << "maxIntensity " << maxIntensity << " minIntensity " << minIntensity << std::endl;  
	//	Fl_Text_Buffer *buff3 = new Fl_Text_Buffer();
	//	sprintf(texto3,"Min %f, Max %f",minIntensity,maxIntensity);
	//	buff3->text(texto3);
	//	m_MaxMinText->buffer(buff3);
	
	panelMessage->show();
}

void TensorConsole::ComputeScalarValue( unsigned int scalar, const char* nameScalar, unsigned int window )
{
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0){return;}
	
	ComputeScalarsPointer filter = ComputeScalarsType::New();
	filter->SetInput( (*m_VectorTensorData)[dataId].image );
	switch( scalar ){
		case 0:
			filter->SetComputeFA();
			break;
		case 1:
			filter->SetComputeRA();
			break;
		case 2:
			filter->SetComputeEigVal(0);
			break;
		case 3:
			filter->SetComputeEigVal(1);
			break;
		case 4:
			filter->SetComputeEigVal(2);
			break;
		case 5:
			filter->SetComputeMD();
			break;
		case 6:
			filter->SetComputeDC(0);
			break;
		case 7:
			filter->SetComputeDC(1);
			break;
		case 8:
			filter->SetComputeDC(2);
			break;
		case 9:
			filter->SetComputeDC(3);
			break;
		case 10:
			filter->SetComputeDC(4);
			break;
		case 11:
			filter->SetComputeDC(5);
			break;
		case 12:
			filter->SetComputeShapeCoefficients(0);
			break;
		case 13:
			filter->SetComputeShapeCoefficients(1);
			break;
		case 14:
			filter->SetComputeShapeCoefficients(2);
			break;
		default:
			filter->SetComputeFA();
			break;
	}
	try{
		filter->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	char nombre[200];
	sprintf( nombre, nameScalar );
	
	m_VectorData->generateNewData( filter->GetOutput(), nombre );
	ImageViewer[window]->SetImage( filter->GetOutput() );
	filter = NULL;
	return;
}

void TensorConsole::ViewStrainSlice3D() {
	ViewStrainSlice3D(this->tiempo);
}

void TensorConsole::ViewStrainSlice3D(int tiempo)
{
	int dataId = m_strainTensorDataBrowser->value()-1;
	if( dataId<0){return;}

	DeformImageType::Pointer image = (*m_VectorSTData)[dataId].deform_image;

	ExtractFilterType::Pointer extract = ExtractFilterType::New();
	extract->SetInput( image );

	DeformImageType::RegionType inputRegion = image->GetLargestPossibleRegion();

	DeformImageType::SizeType inputSize = inputRegion.GetSize();
	inputSize[3] = 0;

	DeformImageType::IndexType inputIndex = inputRegion.GetIndex();
	inputIndex[3] = tiempo;

	DeformImageType::RegionType desiredRegion;
	desiredRegion.SetSize( inputSize );
	desiredRegion.SetIndex( inputIndex );

	extract->SetExtractionRegion( desiredRegion );

	ComputeDeformFilterPointer filter = ComputeDeformFilterType::New();
	filter->SetInput( extract->GetOutput() );

	try{
		filter->Update();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
  	m_VTKexporter->SetInput(filter->GetOutput());

  	vtkImageImport* vtkImporter = vtkImageImport::New();
  	ConnectPipelines(m_VTKexporter, vtkImporter);

  	try{
		  vtkImporter->GetOutput()->Update();
  	}
  	catch(itk::ExceptionObject &excp){
	  	fl_alert( excp.GetDescription() );
	  	return;
  	}

  	this->m_activeStrainImage = vtkImporter->GetOutput();

  	this->m_activeStrainImage->Update();

  	ImageViewerStrain3D->renderSlice( this->m_activeStrainImage, 2, ImageViewerStrain3D->planeWidgetZ->GetSliceIndex() );
  
  	Fl::check();
  	ImageViewerStrain3D->redraw();
  	Fl::check();

	filter = NULL;

	return;
}

void TensorConsole::CreateStreamLineFloat(void) {
	
	std::cout << "Creating streamline " << std::endl;
	int numberOfPoints = 0;
	NodeType node;
	vtkPolyData* streamlineList=vtkPolyData::New();
	streamlineList->Allocate();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	
	//// Recorro las semillas: 
	if (m_seeds) {
		for (unsigned int i=0;i<m_seeds->size();i++) {
			vtkIdList* streamline=vtkIdList::New();
			streamline->Allocate(0);
			std::vector< InputPointType > VectorElementFloat = m_points_out_float[i];
			// Recorro los puntos de cada fibra 
			for ( unsigned j = 0; j < VectorElementFloat.size(); j++ ) {
				points->InsertNextPoint(VectorElementFloat[j][0],VectorElementFloat[j][1],VectorElementFloat[j][2]);
				streamline->InsertNextId(numberOfPoints);
				numberOfPoints++;
				//std::cout << "Semilla " << i << " Punto  " << VectorElementFloat[j][0] << " " << VectorElementFloat[j][1] << " " << VectorElementFloat[j][2] << " numberOfPoints " << numberOfPoints << std::endl;
			} // end for j
			streamlineList->InsertNextCell(VTK_POLY_LINE, streamline);
			streamline->Delete();
		} // end for i
		//std::cout << "Creada streamline de  " << numberOfPoints << " points " << std::endl;
	} 
	///////// 
	streamlineList->SetPoints(points);
	
	vtkActor *Actor = vtkActor::New();
	DataModelElementType model_data;
	model_data.Id     = m_VectorModelData->size();
	model_data.nombre = "New streamline";
	model_data.data   = streamlineList;
	model_data.actor  = Actor;
	m_VectorModelData->push_back(model_data);
	
	/**
	 m_modeldataBrowser->add("Nuevo streamline");
	 m_modeldataBrowser->select(m_modeldataBrowser->size());
	 m_modeldataBrowser->redraw();
	 */
	
	double color[3] = {m_color_r->value(),m_color_g->value(),m_color_b->value()};
	ImageViewer3D->renderTracts( streamlineList, Actor, color, m_radius->value() );
	
	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();
}

void TensorConsole::CreateStreamLineStraight(void) {
	
	std::cout << "Creating streamline " << std::endl;
	int numberOfPoints = 0;
	NodeType node;
	vtkPolyData* streamlineList=vtkPolyData::New();
	streamlineList->Allocate();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	
	//// Recorro las semillas: 
	//std::cout << "Num Semillas: " << m_seeds->size() << std::endl;
	if (m_seeds) {
		for (unsigned int i=0;i<m_seeds->size();i++) {
			vtkIdList* streamline=vtkIdList::New();
			streamline->Allocate(0);
			std::vector< FloatImageType::IndexType> VectorElement = m_points_out[i];
			// Recorro los puntos de cada fibra 
			for ( unsigned j = 0; j < VectorElement.size(); j=j+VectorElement.size()-1 ) {
				points->InsertNextPoint(VectorElement[j][0],VectorElement[j][1],VectorElement[j][2]);
				streamline->InsertNextId(numberOfPoints);
				numberOfPoints++;
			} // end for j
			streamlineList->InsertNextCell(VTK_POLY_LINE, streamline);
			streamline->Delete();
		} // end for i
		//std::cout << "Creada streamline de  " << numberOfPoints << " points " << std::endl;
	} 
	///////// 
	streamlineList->SetPoints(points);
	
	vtkActor *Actor = vtkActor::New();
	DataModelElementType model_data;
	model_data.Id     = m_VectorModelData->size();
	model_data.nombre = "New streamline";
	model_data.data   = streamlineList;
	model_data.actor  = Actor;
	m_VectorModelData->push_back(model_data);
	
	/**
	 m_modeldataBrowser->add("Nuevo streamline");
	 m_modeldataBrowser->select(m_modeldataBrowser->size());
	 m_modeldataBrowser->redraw();
	 */
	
	double color[3] = {0,0,1};
	ImageViewer3D->renderTracts(streamlineList, Actor, color, m_radius->value());
	
	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();
}

void TensorConsole::CreateStreamLine(void) {
	
	std::cout << "Creating streamline " << std::endl;
	int numberOfPoints = 0;
	NodeType node;
	FloatImageType::IndexType index;
	vtkPolyData* streamlineList=vtkPolyData::New();
	streamlineList->Allocate();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	
	vtkCardinalSpline *aSplineX = vtkCardinalSpline::New();
	vtkCardinalSpline *aSplineY = vtkCardinalSpline::New();
	vtkCardinalSpline *aSplineZ = vtkCardinalSpline::New();
	
	//// Recorro las semillas: 
	//std::cout << "Num Semillas: " << m_seeds->size() << std::endl;
	if (m_seeds) {
		int next = 0;
		for (unsigned int i=0;i<m_seeds->size();i++) {
			vtkIdList* streamline=vtkIdList::New();
			streamline->Allocate(0);
			std::vector< FloatImageType::IndexType> VectorElement = m_points_out[i];
			// Recorro los puntos de cada fibra 
			aSplineX->RemoveAllPoints();
			aSplineY->RemoveAllPoints();
			aSplineZ->RemoveAllPoints();
			
			//std::cout << "VectorElement.size: " <<  VectorElement.size() << std::endl;
			int NumPointsIn = 0;
			for ( unsigned j = 0; j < VectorElement.size(); j=j+3 ) {
				index = VectorElement[j];
				//points->InsertNextPoint(index[0],index[1],index[2]);
				//streamline->InsertNextId(numberOfPoints);
				//numberOfPoints++;
				aSplineX->AddPoint(NumPointsIn,index[0]);
				aSplineY->AddPoint(NumPointsIn,index[1]);
				aSplineZ->AddPoint(NumPointsIn,index[2]);
				//std::cout << j << " " <<  index[0] << " " << index[1] << " " << index[2] << std::endl;
				NumPointsIn++;
				//std::cout << "Semilla " << i << " Punto  " << index << " numberOfPoints " << numberOfPoints << std::endl;
			} // end for j
			
			int NumPointsOut = VectorElement.size()-10;
			
			for (int k= 0; k< NumPointsOut; k++) {
				float t = k * ( (float)NumPointsIn - 1.0 ) / ( (float)NumPointsOut - 1.0 );
				points->InsertPoint(next,aSplineX->Evaluate(t),aSplineY->Evaluate(t),aSplineZ->Evaluate(t));
				next++;
				//std::cout << "t " << t << "     " << aSplineX->Evaluate(t) << " " << aSplineY->Evaluate(t) << " " << aSplineZ->Evaluate(t) << std::endl;
				streamline->InsertNextId(numberOfPoints);
				numberOfPoints++;
			}
			
			streamlineList->InsertNextCell(VTK_POLY_LINE, streamline);
			streamline->Delete();
		} // end for i
		//std::cout << "Creada streamline de  " << numberOfPoints << " points " << std::endl;
	} 
	///////// 
	streamlineList->SetPoints(points);
	
	vtkActor *Actor = vtkActor::New();
	DataModelElementType model_data;
	model_data.Id     = m_VectorModelData->size();
	model_data.nombre = "New streamline";
	model_data.data   = streamlineList;
	model_data.actor  = Actor;
	m_VectorModelData->push_back(model_data);
	
	/*
	 m_modeldataBrowser->add("Nuevo streamline");
	 m_modeldataBrowser->select(m_modeldataBrowser->size());
	 m_modeldataBrowser->redraw();
	 */
	
	double color[3] = {1,0,1};
	ImageViewer3D->renderTracts(streamlineList, Actor, color,m_radius->value());
	
	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();
}

void TensorConsole::CreateStreamLineNew(void) {
	
	std::cout << "Doing createStreamLineNew " << std::endl;
	//std::list<ClickPoint> clicked_points;
	//clicked_points = ImageViewer[m_activeinput]->getClickedPointsStored();
	InitializePoints( );
	std::list< ClickPoint >::const_iterator point = clicked_points.begin();
	m_seeds->Initialize();
	int count = 0;
	while( point != clicked_points.end() ) {
		m_seedPosition[0] = point->x;
		m_seedPosition[1] = point->y;
		m_seedPosition[2] = point->z;
		m_node.SetValue( point->value );
		m_node.SetIndex( m_seedPosition );
		m_seeds->InsertElement(count,m_node);
		count++;
		++point;
	}
	
	vtkPolyData *seeds = vtkPolyData::New(); 
	seeds->Allocate();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	FloatImageType::IndexType index;
	
	if (m_seeds) {
		NodeContainer::ConstIterator pointsIter = m_seeds->Begin();
		NodeContainer::ConstIterator pointsEnd = m_seeds->End();
		for ( ; pointsIter != pointsEnd; ++pointsIter ) {
			NodeType node = pointsIter.Value();
			index = node.GetIndex();
			std::cout << "Seed index " << index << std::endl;
			points->InsertNextPoint(index[0],index[1],index[2]);
		} // end for
		//std::cout << "Creada streamline de  " << numberOfPoints << " points " << std::endl;
	} 
	seeds->SetPoints(points);
	
	typedef itk::GradientImageFilter<InputImageType, float, float> GradientImageFilter;
	GradientImageFilter::Pointer gradient = GradientImageFilter::New();
	int Id = m_dataBrowser->value()-1;
	gradient->SetInput((*m_VectorData)[Id].image);
	gradient->Update();
	
	vtkImageData* data;  
	vtkImageImport* vtkImporter = vtkImageImport::New();
	typedef itk::CovariantVector<float, 3> Vector3DPixelType;
	typedef itk::Image< Vector3DPixelType, 3 >    Vector3DImageType;
	typedef itk::VTKImageExport<Vector3DImageType > VTKexporterVectorType;
	VTKexporterVectorType::Pointer VTKexporterVector = VTKexporterVectorType::New();
	VTKexporterVector->SetInput(gradient->GetOutput());
	//m_VTKexporter->SetInput(m_VectorData[m_op1->value()].image);
	ConnectPipelines(VTKexporterVector, vtkImporter);
	vtkImporter->GetOutput()->Update();
	data = vtkImporter->GetOutput();

	vtkDataArray *da;
	da = data->GetPointData()->GetScalars();
	data->GetPointData()->SetVectors(da);

	vtkRungeKutta4 *rk4 = vtkRungeKutta4::New();
	vtkStreamLine *streamer = vtkStreamLine::New();
	streamer->SetInput(data);
	streamer->SetSource(seeds);
	streamer->SetMaximumPropagationTime(m_prop_time->value());
	streamer->SetIntegrationStepLength(m_int_step->value());
	streamer->SetStepLength(m_step_length->value());
	streamer->SetNumberOfThreads(1);
	//streamer->SetIntegrationDirectionToForward();
	streamer->SetIntegrationDirectionToBackward();
	streamer->VorticityOn();
	streamer->SetIntegrator(rk4);
	streamer->Update();
	
	vtkActor *Actor = vtkActor::New();
	DataModelElementType model_data;
	model_data.Id     = m_VectorModelData->size();
	model_data.nombre = "New streamline";
	model_data.data   = seeds;
	model_data.actor  = Actor;
	m_VectorModelData->push_back(model_data);
	
	/*
	 m_modeldataBrowser->add("New streamline");
	 m_modeldataBrowser->select(m_modeldataBrowser->size());
	 m_modeldataBrowser->redraw();
	 */
	
	ImageViewer3D->renderTractsNew(streamer, Actor);
}	

void TensorConsole::geodesicPath3D(void ) 
{
	//std::list<ClickPoint> clicked_points;
	//clicked_points = ImageViewer[m_activeinput]->getClickedPointsStored();
	InitializePoints( );
	std::list< ClickPoint >::const_iterator point = clicked_points.begin();
	m_seeds->Initialize();
	int count = 0;
	while( point != clicked_points.end() ) {
		m_seedPosition[0] = point->x;
		m_seedPosition[1] = point->y;
		m_seedPosition[2] = point->z;
		m_node.SetValue( point->value );
		m_node.SetIndex( m_seedPosition );
		m_seeds->InsertElement(count,m_node);
		count++;
		++point;
	}
	//std::cout << "Num semillas: " << count << "\n";
	m_geodesicPath3D->SetSeeds( m_seeds );
	int Id = m_dataBrowser->value()-1;
	m_geodesicPath3D->SetInput( (*m_VectorData)[Id].image );
	m_geodesicPath3D->SetH(m_h->value());
	m_geodesicPath3D->SetMetodo((int)m_metodo->value());
	m_geodesicPath3D->SetMetodoPtoIntermedio((int)m_metodo_pto_intermedio->value());
	m_geodesicPath3D->SetSigma(m_sigma->value());
	
	try{
		m_geodesicPath3D->Update();
		m_geodesicPath3D->Modified();
	}
	catch( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	
	m_points_out.clear();
	m_points_out = m_geodesicPath3D->GetPointsOut();
	m_points_out_float.clear();
	m_points_out_float = m_geodesicPath3D->GetPointsOutFloat();
	
}

/*
 void TensorConsole::EulerTractography(void ) 
 {
 std::cout << "Doing Euler tractography " << std::endl;
 //std::list<ClickPoint> clicked_points;
 //clicked_points = ImageViewer[m_activeinput]->getClickedPointsStored();
 InitializePoints( );
 std::list< ClickPoint >::const_iterator point = clicked_points.begin();
 m_seeds->Initialize();
 int count = 0;
 while( point != clicked_points.end() ) {
 m_seedPosition[0] = point->x;
 m_seedPosition[1] = point->y;
 m_seedPosition[2] = point->z;
 m_node.SetValue( point->value );
 m_node.SetIndex( m_seedPosition );
 m_seeds->InsertElement(count,m_node);
 count++;
 ++point;
 }
 
 vtkPolyData *seeds = vtkPolyData::New(); 
 seeds->Allocate();
 vtkPoints* points = vtkPoints::New();
 points->Allocate(0);
 FloatImageType::IndexType index;
 
 if (m_seeds) {
 NodeContainer::ConstIterator pointsIter = m_seeds->Begin();
 NodeContainer::ConstIterator pointsEnd = m_seeds->End();
 for ( ; pointsIter != pointsEnd; ++pointsIter ) {
 NodeType node = pointsIter.Value();
 index = node.GetIndex();
 std::cout << "Seed index " << index << std::endl;
 points->InsertNextPoint(index[0],index[1],index[2]);
 } // end for
 } 
 seeds->SetPoints(points);
 
 vtkImageData* data;  
 vtkImageImport* vtkImporter = vtkImageImport::New();
 
 typedef itk::ImageRegionIteratorWithIndex< TensorImageType>				TensorIteratorType;
 TensorIteratorType	tensorIt((*m_VectorTensorData)[0].image, (*m_VectorTensorData)[0].image->GetRequestedRegion());
 
 typedef itk::Vector<double,3>                         EigenVectorType;
 typedef itk::Image<EigenVectorType,3>                  EigenvectorImageType;
 
 EigenvectorImageType::Pointer eigenvector=EigenvectorImageType::New();
 eigenvector->SetRegions((*m_VectorTensorData)[0].image->GetLargestPossibleRegion());
 eigenvector->Allocate();
 typedef itk::ImageRegionIteratorWithIndex< EigenvectorImageType>						VectorIteratorType;
 VectorIteratorType	vectorIt(eigenvector, eigenvector->GetRequestedRegion());
 
 typedef itk::DTITensor<float>::EigenValuesArrayType                     EigenValuesArrayType;
 EigenValuesArrayType eigval;
 double	eigvec[3];
 for (tensorIt.GoToBegin(),vectorIt.GoToBegin();!tensorIt.IsAtEnd();++tensorIt,++vectorIt) {
 
 tensorIt.Get().ComputeEigenValues( eigval );
 tensorIt.Get().ComputeEigenVector( eigval[0], eigvec);
 vectorIt.Set(eigvec);
 }
 
 typedef itk::Vector<double, 3> Vector3Type;
 typedef itk::Image< Vector3Type, 3 >    Vector3ImageType;
 typedef itk::VTKImageExport< Vector3ImageType> VTKexportVector3Type;
 VTKexportVector3Type::Pointer m_VTKexporter3Vector = VTKexportVector3Type::New(); 
 
 m_VTKexporter3Vector->SetInput(eigenvector);
 ConnectPipelines(m_VTKexporter3Vector, vtkImporter);
 vtkImporter->GetOutput()->Update();
 data = vtkImporter->GetOutput();
 
 vtkDataArray *da;
 da = data->GetPointData()->GetScalars();
 data->GetPointData()->SetVectors(da);
 
 vtkEuler *rk4 = vtkEuler::New();
 vtkStreamLine *streamer = vtkStreamLine::New();
 streamer->SetInput(data);
 streamer->SetSource(seeds);
 streamer->SetMaximumPropagationTime(m_prop_time->value());
 streamer->SetIntegrationStepLength(m_int_step->value());
 streamer->SetStepLength(m_step_length->value());
 streamer->SetNumberOfThreads(1);
 //streamer->SetIntegrationDirectionToForward();
 streamer->SetIntegrationDirectionToBackward();
 streamer->VorticityOn();
 streamer->SetIntegrator(rk4);
 streamer->Update();
 
 vtkActor *Actor = vtkActor::New();
 DataModelElementType model_data;
 model_data.Id = (*m_VectorModelData).size();
 model_data.nombre = "New streamline";
 model_data.data = seeds;
 model_data.actor = Actor;
 (*m_VectorModelData).push_back(model_data);
 
 m_modeldataBrowser->add("New streamline");
 m_modeldataBrowser->select(m_modeldataBrowser->size());
 m_modeldataBrowser->redraw();
 
 ImageViewer3D->renderTractsNew(streamer, Actor);
 }*/

void TensorConsole::InitializePoints(  ) 
{
	
	clicked_points.clear(); 
	
	UCharIteratorType It(m_imageOverlay ,m_imageOverlay->GetRequestedRegion());
	for (It.GoToBegin();!It.IsAtEnd();++It) {
		if (It.Get() != 0) {
			//std::cout << "Index: " << It.GetIndex() << " Valor: " << (unsigned char)It.Get() << std::endl;
		    ClickPoint clickedPoint(It.GetIndex()[0], It.GetIndex()[1], It.GetIndex()[2], 0, It.Get());
			clicked_points.push_front( clickedPoint );
		}
	}
	
}

void TensorConsole::InitializePoints( unsigned int index ) 
{
	
	clicked_points.clear(); 
	
	UCharIteratorType It(m_imageOverlay ,m_imageOverlay->GetRequestedRegion());
	for (It.GoToBegin();!It.IsAtEnd();++It) {
		if (It.Get() == index) {
		    ClickPoint clickedPoint(It.GetIndex()[0], It.GetIndex()[1], It.GetIndex()[2], 0, It.Get());
			clicked_points.push_front( clickedPoint );
		}
	}
	
}

void TensorConsole::RetainNearbyFibers(vtkPolyData* streamLineList, float factor, int labelROI) {
	vtkPolyData* NewstreamLineList = vtkPolyData::New();
	NewstreamLineList->Allocate();
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points_fibra = vtkPoints::New();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	points_fibra->Allocate(0);
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	std::vector<float> vector_distancias;
	float mean_size = 0;
	double punto[3],punto_ant[3];
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points_fibra = polyline->GetPoints();		
        mean_size += points_fibra->GetNumberOfPoints();
	}
	mean_size = mean_size/streamLineList->GetNumberOfCells();
	//std::cout << "Mean number of points: " << mean_size << std::endl;
	int numberOfPoints = 0;
	// para cada fibra
	float dist = 999999;
	for (int fibra=0;fibra<streamLineList->GetNumberOfCells();fibra++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(fibra);
		points_fibra = polyline->GetPoints();
		vtkIdList* streamline=vtkIdList::New();
		streamline->Allocate(0);
		int numero_de_puntos = points_fibra->GetNumberOfPoints();
		vector_distancias.resize(numero_de_puntos);
		//forward computation:
		for (int k=0;k< numero_de_puntos ;k++) {
			points_fibra->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			if (m_imageOverlay->GetPixel(index) == labelROI) {
				dist = 0;
			} else {
			  if (k>0) { 
			    points_fibra->GetPoint(k-1,punto_ant);
			    dist += sqrt((punto_ant[0] - punto[0])*(punto_ant[0] - punto[0]) + 
						(punto_ant[1] - punto[1])*(punto_ant[1] - punto[1]) +
						(punto_ant[2] - punto[2])*(punto_ant[2] - punto[2])); 
			  }
			  vector_distancias[k] = dist;
			}	
		}
		//reverse computation
		dist = 999999;
		int nuevo_numero_de_puntos = 0;
		for (int k=numero_de_puntos-1;k>=0;k--) {
			points_fibra->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			if (m_imageOverlay->GetPixel(index) == labelROI) {
				dist = 0;
			} else {
			  if (k< numero_de_puntos-1) { 
			    points_fibra->GetPoint(k+1,punto_ant);
				dist += sqrt((punto_ant[0] - punto[0])*(punto_ant[0] - punto[0]) + 
						(punto_ant[1] - punto[1])*(punto_ant[1] - punto[1]) +
						(punto_ant[2] - punto[2])*(punto_ant[2] - punto[2])); 
			  }
			}
			if (dist < vector_distancias[k]) { 
			  vector_distancias[k] = dist;
			}
			if (vector_distancias[k] < mean_size*factor) { 
			  points->InsertNextPoint(punto);
			  streamline->InsertNextId(numberOfPoints);
			  nuevo_numero_de_puntos++;
			  numberOfPoints++;
			}
		}
		//std::cout << "mean size*factor " << mean_size*factor << std::endl; 
		//std::cout << "vector distancias " << std::endl; 
		//for (int i=0;i<numero_de_puntos;i++) std::cout << vector_distancias[i] << " ";
	    //std::cout << "numero de puntos de fibra inicial: " << numero_de_puntos << " final: " << nuevo_numero_de_puntos << std::endl;

		NewstreamLineList->InsertNextCell(VTK_POLY_LINE,streamline); 
		streamline->Delete();
	}
	NewstreamLineList->SetPoints(points);
	// Ahora cambio uno por otro: 
	streamLineList->DeepCopy(NewstreamLineList);	
}

// Ordena las fibras en el eje "y" (de anterior a posterior) en un  plano z (axial) determinado por el usuario
void TensorConsole::OrderFibersCoronal(vtkPolyData* streamLineList, int z) {
	vtkPolyData* NewstreamLineList = vtkPolyData::New();
	NewstreamLineList->Allocate();
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points_fibra = vtkPoints::New();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	points_fibra->Allocate(0);
	TensorImageType::IndexType index;
	double punto[3];
	std::list<int> L,L2;
	std::list< int >::iterator Lit,Lit2;
	double Z_primero,Z_ultimo;
	int numberOfPoints = 0;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
		
	// recorremos todas las fibras para ordenarlas en direccion A-P
	for (int fibra=0;fibra<streamLineList->GetNumberOfCells();fibra++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(fibra);
		points_fibra = polyline->GetPoints();
		//std::cout << "fibra " << fibra << " L.size() " << L.size() << std::endl; 
		int numero_de_puntos = points_fibra->GetNumberOfPoints();
		// Para los ptos de la fibra
		for (int k=0;k< numero_de_puntos ;k++) {
			points_fibra->GetPoint(k,punto);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			if ( index[2] == z) {
				if (L.size() == 0) {
				  L.push_back(fibra);
				  L2.push_back(index[1]);
				} else {
				  int flag = 0;
				  for (Lit= L.begin(), Lit2 = L2.begin();Lit != L.end();++Lit,++Lit2) {
					if (*Lit2 > index[1]) {
					  L.insert(Lit, fibra);
					  L2.insert(Lit2,index[1]);
					  flag = 1;
					  break;	
					}
				  }
				  if (flag == 0) {
					  L.push_back(fibra);
					  L2.push_back(index[1]);
				  }
				}
				break; 
			}
		}
	}

	
	// Recorremos todas las fibras de forma ordenada en la direccion y (perpendicular al plano coronal) 
	// para almacenarlas en orden A-P
	// Tambien forzamos que los puntos de cada fibra esten ordenados en el mismo sentido en el eje z, en direccion I-S
	for (Lit= L.begin();Lit != L.end();++Lit) {
		//std::cout << "Añadiendo fibra " << *Lit << std::endl; 
		polyline = (vtkPolyLine*)streamLineList->GetCell(*Lit);
		points_fibra = polyline->GetPoints();
		int numero_de_puntos = points_fibra->GetNumberOfPoints();
		vtkIdList* streamline=vtkIdList::New();
		streamline->Allocate(0);
		// Comprobamos que todos los puntos vayan en el mismo sentido: 
		points_fibra->GetPoint(0,punto);
		Z_primero = punto[2];
		points_fibra->GetPoint(numero_de_puntos-1,punto);
		Z_ultimo = punto[2];
		// Para los ptos de la fibra
		if ((Z_primero - Z_ultimo) < 0) {
		  //std::cout << "hay que ordenar al reves, Primero: " << Z_primero << " ultimo " << Z_ultimo << std::endl;
			for (int k=numero_de_puntos-1;k>=0 ;k--) {
				points_fibra->GetPoint(k,punto);
				points->InsertNextPoint(punto);
				streamline->InsertNextId(numberOfPoints);
				numberOfPoints++;
				newScalars->InsertNextValue(*Lit);
		  }
		} else {
		  for (int k=0;k< numero_de_puntos ;k++) {
			points_fibra->GetPoint(k,punto);
			points->InsertNextPoint(punto);
			streamline->InsertNextId(numberOfPoints);
			numberOfPoints++;
			newScalars->InsertNextValue(*Lit);
		  }
		}
		NewstreamLineList->InsertNextCell(VTK_POLY_LINE,streamline); 
		streamline->Delete();
	}
	
	
	NewstreamLineList->SetPoints(points);
	if (newScalars->GetNumberOfTuples() > 0) {
		NewstreamLineList->GetPointData()->SetScalars(newScalars);
		int dataId = m_modeldataBrowser->value()-1;
		if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
		(*m_VectorModelData)[dataId].actor->GetMapper()->SetScalarRange(0,L.size());
	}
	
	// Ahora cambio uno por otro: 
	streamLineList->DeepCopy(NewstreamLineList);	
}

// Ordena las fibras en el eje "y" (de anterior a posterior) fijando un plano y (sagital) determinado por el usuario
void TensorConsole::OrderFibersCoronal2(vtkPolyData* streamLineList, int y) {
	vtkPolyData* NewstreamLineList = vtkPolyData::New();
	NewstreamLineList->Allocate();
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points_fibra = vtkPoints::New();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	points_fibra->Allocate(0);
	TensorImageType::IndexType index;
	double punto[3];
	std::list<int> L,L2;
	std::list< int >::iterator Lit,Lit2;
	double X_primero,X_ultimo;
	int numberOfPoints = 0;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	// recorremos todas las fibras para ordenarlas en direccion L-R
	for (int fibra=0;fibra<streamLineList->GetNumberOfCells();fibra++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(fibra);
		points_fibra = polyline->GetPoints();
		//std::cout << "fibra " << fibra << " L.size() " << L.size() << std::endl; 
		int numero_de_puntos = points_fibra->GetNumberOfPoints();
		// Para los ptos de la fibra
		for (int k=0;k< numero_de_puntos ;k++) {
			points_fibra->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);			
			if ( index[0] == y || index[0] == y+1 || index[0] == y-1) {
				if (L.size() == 0) {
					L.push_back(fibra);
					L2.push_back(index[1]);
				} else {
					int flag = 0;
					for (Lit= L.begin(), Lit2 = L2.begin();Lit != L.end();++Lit,++Lit2) {
						if (*Lit2 > index[1]) {
							L.insert(Lit, fibra);
							L2.insert(Lit2,index[1]);
							flag = 1;
							break;	
						}
					}
					if (flag == 0) {
						L.push_back(fibra);
						L2.push_back(index[1]);
					}
				}
				break; 
			}
		}
	}
		
	// Recorremos todas las fibras de forma ordenada en la direccion y (perpendicular al plano coronal) 
	// para almacenarlas en orden A-P
	// Tambien forzamos que los puntos de cada fibra esten ordenados en el mismo sentido en el eje z, en direccion I-S
	for (Lit= L.begin();Lit != L.end();++Lit) {
		//std::cout << "Añadiendo fibra " << *Lit << std::endl; 
		polyline = (vtkPolyLine*)streamLineList->GetCell(*Lit);
		points_fibra = polyline->GetPoints();
		int numero_de_puntos = points_fibra->GetNumberOfPoints();
		vtkIdList* streamline=vtkIdList::New();
		streamline->Allocate(0);
		// Comprobamos que todos los puntos vayan en el mismo sentido: 
		points_fibra->GetPoint(0,punto);
		X_primero = punto[0];
		points_fibra->GetPoint(numero_de_puntos-1,punto);
		X_ultimo = punto[0];
		// Para los ptos de la fibra
		if ((X_primero - X_ultimo) < 0) {
			//std::cout << "hay que ordenar al reves, Primero: " << Z_primero << " ultimo " << Z_ultimo << std::endl;
			for (int k=numero_de_puntos-1;k>=0 ;k--) {
				points_fibra->GetPoint(k,punto);
				points->InsertNextPoint(punto);
				streamline->InsertNextId(numberOfPoints);
				numberOfPoints++;
				newScalars->InsertNextValue(*Lit);
			}
		} else {
			for (int k=0;k< numero_de_puntos ;k++) {
				points_fibra->GetPoint(k,punto);
				points->InsertNextPoint(punto);
				streamline->InsertNextId(numberOfPoints);
				numberOfPoints++;
				newScalars->InsertNextValue(*Lit);
			}
		}
		NewstreamLineList->InsertNextCell(VTK_POLY_LINE,streamline); 
		streamline->Delete();
	}
	
	
	NewstreamLineList->SetPoints(points);
	if (newScalars->GetNumberOfTuples() > 0) {
		NewstreamLineList->GetPointData()->SetScalars(newScalars);
		int dataId = m_modeldataBrowser->value()-1;
		if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
		(*m_VectorModelData)[dataId].actor->GetMapper()->SetScalarRange(0,L.size());
	}
	
	// Ahora cambio uno por otro: 
	streamLineList->DeepCopy(NewstreamLineList);	
}

void TensorConsole::CutFibers(vtkPolyData* streamLineList, int orientation, int value, int direccion) {
	vtkPolyData* NewstreamLineList = vtkPolyData::New();
	NewstreamLineList->Allocate();
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points_fibra = vtkPoints::New();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	points_fibra->Allocate(0);
	TensorImageType::IndexType index;
	
	if (orientation > 2 || orientation < 0) return;
	
	double punto[3];
	int numberOfPoints = 0;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	// para cada fibra
	for (int fibra=0;fibra<streamLineList->GetNumberOfCells();fibra++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(fibra);
		points_fibra = polyline->GetPoints();
		int numero_de_puntos = points_fibra->GetNumberOfPoints();
		vtkIdList* streamline=vtkIdList::New();
		streamline->Allocate(0);
		// Para los ptos de la fibra
		for (int k=0;k< numero_de_puntos ;k++) {
			points_fibra->GetPoint(k,punto);
			index[orientation] = round((punto[orientation]-origen[orientation])/spacing[orientation]);
			if (direccion == 0) {
			  if ( index[orientation] >= value) {
				points->InsertNextPoint(punto);
				streamline->InsertNextId(numberOfPoints);
				numberOfPoints++;
			  }
			}
			if (direccion == 1) {
				if ( index[orientation] <= value) {
					points->InsertNextPoint(punto);
					streamline->InsertNextId(numberOfPoints);
					numberOfPoints++;
				}
			}
		}
		NewstreamLineList->InsertNextCell(VTK_POLY_LINE,streamline); 
		streamline->Delete();
	}
	NewstreamLineList->SetPoints(points);
	// Ahora cambio uno por otro: 
	streamLineList->DeepCopy(NewstreamLineList);	
	
}

void TensorConsole::RemoveLargeFibers(vtkPolyData* streamLineList, float factor) {
	vtkPolyData* NewstreamLineList = vtkPolyData::New();
	NewstreamLineList->Allocate();
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points_fibra = vtkPoints::New();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	points_fibra->Allocate(0);
	float mean_size = 0;
	double punto[3];
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points_fibra = polyline->GetPoints();		
        mean_size += points_fibra->GetNumberOfPoints();
	}
	mean_size = mean_size/streamLineList->GetNumberOfCells();
	//std::cout << "Mean number of points: " << mean_size << std::endl;
	int numberOfPoints = 0;
	// para cada fibra
	for (int fibra=0;fibra<streamLineList->GetNumberOfCells();fibra++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(fibra);
		points_fibra = polyline->GetPoints();
		int numero_de_puntos = points_fibra->GetNumberOfPoints();
		if ((float)numero_de_puntos < factor*mean_size ) {
		  vtkIdList* streamline=vtkIdList::New();
		  streamline->Allocate(0);
		  // Para los ptos de la fibra
		  for (int k=0;k< numero_de_puntos ;k++) {
			points_fibra->GetPoint(k,punto);
			points->InsertNextPoint(punto);
			streamline->InsertNextId(numberOfPoints);
			numberOfPoints++;
		  }
		  NewstreamLineList->InsertNextCell(VTK_POLY_LINE,streamline); 
		  streamline->Delete();
		}
	}
	NewstreamLineList->SetPoints(points);
	// Ahora cambio uno por otro: 
	streamLineList->DeepCopy(NewstreamLineList);	
}

void TensorConsole::ResizeLargeFibers(vtkPolyData* streamLineList, float factor) {
	vtkPolyData* NewstreamLineList = vtkPolyData::New();
	NewstreamLineList->Allocate();
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points_fibra = vtkPoints::New();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	points_fibra->Allocate(0);
	float mean_size = 0;
	double punto[3];
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points_fibra = polyline->GetPoints();		
        mean_size += points_fibra->GetNumberOfPoints();
	}
	mean_size = mean_size/streamLineList->GetNumberOfCells();
	//std::cout << "Mean number of points: " << mean_size << std::endl;
	int numberOfPoints = 0;
	// para cada fibra
	for (int fibra=0;fibra<streamLineList->GetNumberOfCells();fibra++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(fibra);
		points_fibra = polyline->GetPoints();
		vtkIdList* streamline=vtkIdList::New();
		streamline->Allocate(0);
		int numero_de_puntos = points_fibra->GetNumberOfPoints();
		// Para los ptos de la fibra
		for (int k=0;k< numero_de_puntos ;k++) {
		  if ((float)k < factor*mean_size ) {
			 points_fibra->GetPoint(k,punto);
			 points->InsertNextPoint(punto);
			 streamline->InsertNextId(numberOfPoints);
			 numberOfPoints++;
		  }
		}
		NewstreamLineList->InsertNextCell(VTK_POLY_LINE,streamline); 
		streamline->Delete();
	}
	NewstreamLineList->SetPoints(points);
	// Ahora cambio uno por otro: 
	streamLineList->DeepCopy(NewstreamLineList);	
}

void TensorConsole::ComputeFiberMeasures(bool save_stats_file) {
  int dataId = m_modeldataBrowser->value()-1;
  if (dataId < 0) return;
  for (int i=0;i<19;i++) m_buff[i]->text("");
  if (m_Fiber_measure[0]->value() == 1) this->ComputeFibersFA((*m_VectorModelData)[dataId].data, false, save_stats_file);
  if (m_Fiber_measure[1]->value() == 1) this->ComputeFibersRA((*m_VectorModelData)[dataId].data, save_stats_file);
  if (m_Fiber_measure[2]->value() == 1) this->ComputeFibersMD((*m_VectorModelData)[dataId].data, save_stats_file);
  if (m_Fiber_measure[3]->value() == 1) this->ComputeFibersSize((*m_VectorModelData)[dataId].data, save_stats_file);
  if (m_Fiber_measure[7]->value() == 1 ) this->ComputeFibersEigen((*m_VectorModelData)[dataId].data, save_stats_file);
  if (m_Fiber_measure[8]->value() == 1 ) this->ComputeFibersShapeCoefficients((*m_VectorModelData)[dataId].data, save_stats_file);
  if (m_Fiber_measure[9]->value() == 1 ) this->ComputeFibersTensorComponents((*m_VectorModelData)[dataId].data, save_stats_file);
  this->ShowMeasures();
}
	
void TensorConsole::ComputeFibersSize(vtkPolyData* streamLineList,  bool save_stats_file) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float mean_size = 0, min_size = 10000000, max_size = 0;
	std::vector<float> vec_sizes,vec_curvatures;
	double punto[3],punto_ant[3];
	std::ofstream Sfile;
	
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		float dist = 0;
		
		
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			if (k > 0 && k<points->GetNumberOfPoints()) {
				points->GetPoint(k,punto);
				points->GetPoint(k-1,punto_ant);
				dist += sqrt((punto_ant[0] - punto[0])*(punto_ant[0] - punto[0]) + 
			  		         (punto_ant[1] - punto[1])*(punto_ant[1] - punto[1]) +
							 (punto_ant[2] - punto[2])*(punto_ant[2] - punto[2])); 
			}
		}
		vec_sizes.push_back(dist);
		if (dist > max_size) max_size = dist;
		if (dist < min_size) min_size = dist;
		mean_size += dist;
		//std::cout << "Fiber " << i << " , Size = " << dist << std::endl;
	}
	mean_size = mean_size/streamLineList->GetNumberOfCells();
	
	//std::cout << "-- Num of fibers computed: " << streamLineList->GetNumberOfCells()  << std::endl;
//	std::cout << "Mean size of fibers: " << mean_size << std::endl;
//	std::cout << "Max size of fibers: " << max_size << std::endl;
//	std::cout << "Min size of fibers: " << min_size << std::endl;
	char st[200]; 
	sprintf(st,"Size \nMean %f\nMax %f\nMin %f",mean_size,max_size,min_size);
	m_buff[3]->text(st);
	if (save_stats_file) {
		Sfile.open(m_filename, ios::app);
		Sfile<< st << "\n------------" << std::endl;
		Sfile.close();
	}
}

void TensorConsole::ComputeFibersAngle(vtkPolyData* streamLineList) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float mod_vn,mod_vn_ant,c_aux,curvature,mean_curvature, mean_global_curvature = 0, min_curvature = 10000, max_curvature = 0;
	std::vector<float> vec_curvatures;
	itk::Vector<float,3> vn,vn_ant;
	double punto[3],punto_ant[3],punto_post[3];
	
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		mean_curvature = 0; 
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			if (k > 0 && k<points->GetNumberOfPoints()) {
				points->GetPoint(k,punto);
				points->GetPoint(k-1,punto_ant);
				
				vn_ant[0] = punto[0] - punto_ant[0];
				vn_ant[1] = punto[1] - punto_ant[1];
				vn_ant[2] = punto[2] - punto_ant[2];
			}
			if (k> 0 && k<points->GetNumberOfPoints()-1) {
				points->GetPoint(k+1,punto_post);
			    vn[0] = punto_post[0] - punto[0];
				vn[1] = punto_post[1] - punto[1];
				vn[2] = punto_post[2] - punto[2];
			    c_aux=fabs(vn*vn_ant);
			    mod_vn=sqrt(vn*vn);
			    mod_vn_ant=sqrt(vn_ant*vn_ant);				
			    c_aux=c_aux/(mod_vn*mod_vn_ant);
			    if(c_aux>1) c_aux = 1;
			    curvature=fabs(acos( c_aux ));
				mean_curvature += curvature;
			    //std::cout << "Punto: " << punto[0] << " " << punto[1] << " " << punto[2] << std::endl;
			    //std::cout << "Angle: " << curvature << std::endl;
			}
		}
		if (points->GetNumberOfPoints()> 2) {
			mean_curvature = mean_curvature/(points->GetNumberOfPoints()-2);
		} else {
			mean_curvature = 0;
		}
		mean_global_curvature += mean_curvature;
		vec_curvatures.push_back(mean_curvature);
		if (mean_curvature > max_curvature) max_curvature = mean_curvature;
		if (mean_curvature < min_curvature) min_curvature = mean_curvature;
		//std::cout << "Fiber " << i << " , Mean angle = " << mean_curvature*180/3.1416 << std::endl;
	}
	mean_global_curvature = mean_global_curvature/streamLineList->GetNumberOfCells();
	std::cout << "-- Num of fibers computed: " << streamLineList->GetNumberOfCells()  << std::endl;
	std::cout << "Mean global angle: " << mean_global_curvature*180/3.1416 << std::endl;
	std::cout << "Max angle: " << max_curvature*180/3.1416 << std::endl;
	std::cout << "Min angle: " << min_curvature*180/3.1416 << std::endl;
}

void TensorConsole::ComputeFibersFA(vtkPolyData* streamLineList, bool save_profile, bool save_stats_file) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float mean_FA_tot, mean_FA, max_FA = 0, min_FA = 10, mean_global_FA = 0, total_FA = 0;
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	std::ofstream FAfile;
	std::cout << "Computing FA" << std::endl;
	if (save_profile) {
	  FAfile.open("Profile.txt");
	}
	int z_cut = 10;
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		mean_FA_tot = 0;
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=points->GetNumberOfPoints()-1;k>=0;k--) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			mean_FA = (float)aux->GetPixel(index).GetFractionalAnisotropy();
			// Si se va a grabar el fichero, estamos al principio de la fibra, y la fibra empieza por encima de z_cut rellenamos con ceros al principio 
			if (save_profile && k==points->GetNumberOfPoints()-1 && index[2] > z_cut) {
				//std::cout << "Insertando puntos en la siguente fibra: " << index[2]-z_cut << std::endl;
				points->GetPoint(k-1,punto);
				index[2] = round((punto[2]-origen[2])/spacing[2]);
			    // Momentaneamente comentado para los datos sinteticos:
				//for (int j=0;j<index[2]-z_cut;j++) FAfile<<"0 ";
			}
			if (save_profile) 	FAfile<<mean_FA<<" ";
			mean_FA_tot=mean_FA_tot+mean_FA;
		}
		mean_FA = mean_FA_tot/(points->GetNumberOfPoints());
		if (mean_FA > max_FA) max_FA = mean_FA;
		if (mean_FA < min_FA) min_FA = mean_FA;
		mean_global_FA += mean_FA;
		total_FA += mean_FA_tot;
		//std::cout << "Fiber " << i << " , Mean FA = " << mean_FA << std::endl;
		if (save_profile) FAfile<<std::endl;
	}
    if (save_profile) FAfile.close();
	
	mean_global_FA = mean_global_FA/streamLineList->GetNumberOfCells();
	//std::cout << "-- Num of fibers computed: " << streamLineList->GetNumberOfCells()  << std::endl;
	//std::cout << "Mean FA: "<< mean_global_FA << std::endl;
	//std::cout << "Max FA: " << max_FA<< std::endl;
	//std::cout << "Min FA: " << min_FA << std::endl;
	//std::cout << "Tot FA: " << total_FA << std::endl;	
	char st[200]; 
	sprintf(st,"Mean %f\nMax %f\nMin %f\nTot %f \nI %f",mean_global_FA,max_FA,min_FA,total_FA,total_FA/streamLineList->GetNumberOfCells());
    m_buff[0]->text(st);
	char st2[10];
	sprintf(st2,"%d",streamLineList->GetNumberOfCells());
	m_buffNfibers->text(st2);
	if (save_stats_file) {
		FAfile.open(m_filename);
		FAfile<< "Number of fibers: " << st2 << "\n------------"<< std::endl;
		FAfile<< "FA\n" << st << "\n------------"<< std::endl;
		FAfile.close();
	}
	
}

void TensorConsole::ComputeFibersRA(vtkPolyData* streamLineList, bool save_stats_file) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float mean_RA_tot, mean_RA, max_RA = 0, min_RA = 10, mean_global_RA = 0, total_RA = 0;
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	std::ofstream RAfile;
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		mean_RA_tot = 0;
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=points->GetNumberOfPoints()-1;k>=0;k--) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			mean_RA = (float)aux->GetPixel(index).GetRelativeAnisotropy();
			mean_RA_tot=mean_RA_tot+mean_RA;
		}
		mean_RA = mean_RA_tot/(points->GetNumberOfPoints());
		if (mean_RA > max_RA) max_RA = mean_RA;
		if (mean_RA < min_RA) min_RA = mean_RA;
		mean_global_RA += mean_RA;
		total_RA += mean_RA_tot;
	}
	
	mean_global_RA = mean_global_RA/streamLineList->GetNumberOfCells();
	char st[200]; 
	sprintf(st,"Mean %f\nMax %f\nMin %f\nTot %f",mean_global_RA,max_RA,min_RA,total_RA);
    m_buff[1]->text(st);
	if (save_stats_file) {
		RAfile.open(m_filename, ios::app);
		RAfile<< "RA\n" << st << "\n------------" << std::endl;
		RAfile.close();
	}
}

void TensorConsole::ComputeFibersMD(vtkPolyData* streamLineList, bool save_stats_file) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float mean_MD_tot, mean_MD, max_MD = 0, min_MD = 10, mean_global_MD = 0, total_MD = 0;
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	std::ofstream MDfile;
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		mean_MD_tot = 0;
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=points->GetNumberOfPoints()-1;k>=0;k--) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			mean_MD = (float)aux->GetPixel(index).GetMeanDiffusivity();
			mean_MD_tot=mean_MD_tot+mean_MD;
		}
		mean_MD = mean_MD_tot/(points->GetNumberOfPoints());
		if (mean_MD > max_MD) max_MD = mean_MD;
		if (mean_MD < min_MD) min_MD = mean_MD;
		mean_global_MD += mean_MD;
		total_MD += mean_MD_tot;
	}
	
	mean_global_MD = mean_global_MD/streamLineList->GetNumberOfCells();
	char st[200]; 
	sprintf(st,"Mean %f\nMax %f\nMin %f\nTot %f",mean_global_MD,max_MD,min_MD,total_MD);
    m_buff[2]->text(st);
	if (save_stats_file) {
		MDfile.open(m_filename, ios::app);
		MDfile<< "MD\n" << st << "\n------------" << std::endl;
		MDfile.close();
	}
}

void TensorConsole::ComputeFibersEigen(vtkPolyData* streamLineList, bool save_stats_file) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float  mean_Eig[3], max_Eig[3], min_Eig[3], mean_global_Eig[3], tot_Eig[3];
	for (int j=0;j<3;j++) {
	  mean_global_Eig[j] = 0;
	  min_Eig[j] = 9999;
	  max_Eig[j] = 0;
	  tot_Eig[j] = 0;
	}
	std::ofstream Eigfile;
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=points->GetNumberOfPoints()-1;k>=0;k--) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			TensorPixelType::EigenValuesArrayType eig;
			aux->GetPixel(index).ComputeEigenValues(eig);
			for (int j=0;j<3;j++) mean_Eig[j] = eig[j];
		}
		for (int j=0;j<3;j++) {
		  tot_Eig[j] += mean_Eig[j];
		  mean_Eig[j] = mean_Eig[j]/(points->GetNumberOfPoints());
		  if (mean_Eig[j] > max_Eig[j]) max_Eig[j] = mean_Eig[j];
		  if (mean_Eig[j] < min_Eig[j]) min_Eig[j] = mean_Eig[j];
		  mean_global_Eig[j] += mean_Eig[j];
		}
	}
	
	for (int j=0;j<3;j++) {
	  mean_global_Eig[j] = mean_global_Eig[j]/streamLineList->GetNumberOfCells();
	}
	
	// En realidad min y max son el min y el max de las medias de todas las fibras calculadas
	char st[200]; 
	if (save_stats_file) Eigfile.open(m_filename, ios::app);
	for (int j=0;j<3;j++) {
	  sprintf(st,"Mean %f\nMax %f\nMin %f \nTot %f",mean_global_Eig[j],max_Eig[j],min_Eig[j],tot_Eig[j]);
      m_buff[7+j]->text(st);
	  if (save_stats_file) {
		  Eigfile<< "Eigenvalue "<< j+1 <<"\n" << st << "\n------------" << std::endl;
	  }
	}
	if (save_stats_file) Eigfile.close();
}

void TensorConsole::ComputeFibersShapeCoefficients(vtkPolyData* streamLineList, bool save_stats_file) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	double  cl,cp,cs;
	float  mean_c[3], max_c[3], min_c[3], mean_global_c[3], tot_c[3];
	for (int j=0;j<3;j++) {
	  mean_global_c[j] = 0;
	  min_c[j] = 9999;
	  max_c[j] = 0;
	  tot_c[j] = 0;
	}
	std::ofstream Cfile;
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=points->GetNumberOfPoints()-1;k>=0;k--) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			aux->GetPixel(index).ComputeShapeCoefficients(cl,cp,cs);
			mean_c[0] = cl;
			mean_c[1] = cp;
			mean_c[2] = cs;
		}
		for (int j=0;j<3;j++) {
			tot_c[j] += mean_c[j];
			mean_c[j] = mean_c[j]/(points->GetNumberOfPoints());
			if (mean_c[j] > max_c[j]) max_c[j] = mean_c[j];
			if (mean_c[j] < min_c[j]) min_c[j] = mean_c[j];
			mean_global_c[j] += mean_c[j];
		}
	}
	
	for (int j=0;j<3;j++) {
		mean_global_c[j] = mean_global_c[j]/streamLineList->GetNumberOfCells();
	}
	
	// En realidad min y max son el min y el max de las medias de todas las fibras calculadas
	char st[200]; 
	if (save_stats_file) Cfile.open(m_filename, ios::app);
	char name[3][3] = {"Cl","Cp","Cs"};
	
	for (int j=0;j<3;j++) {
	  sprintf(st,"Mean %f\nMax %f\nMin %f \nTot %f",mean_global_c[j],max_c[j],min_c[j],tot_c[j]);
	  m_buff[10+j]->text(st);
	  if (save_stats_file) {
		  Cfile << name[j] << "\n" << st << "\n------------" << std::endl;
	  }
	}
	if (save_stats_file) Cfile.close();
}

void TensorConsole::ComputeFibersTensorComponents(vtkPolyData* streamLineList, bool save_stats_file) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float  c[6],mean_c[6], max_c[6], min_c[6], mean_global_c[6], tot_c[6];
	for (int i=0;i<6;i++) {
	  mean_global_c[i] = 0;
	  min_c[i] = 9999;
	  max_c[i] = 0;
	  tot_c[i] = 0;
	}
	std::ofstream Dfile;

	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=points->GetNumberOfPoints()-1;k>=0;k--) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			for (int j=0;j<6;j++) {
			  c[j] = aux->GetPixel(index)[j];
			  mean_c[j] = c[j];
			}
		}
		for (int j=0;j<6;j++) {
			tot_c[j] = mean_c[j];
			mean_c[j] = mean_c[j]/(points->GetNumberOfPoints());
			if (mean_c[j] > max_c[j]) max_c[j] = mean_c[j];
			if (mean_c[j] < min_c[j]) min_c[j] = mean_c[j];
			mean_global_c[j] += mean_c[j];
		}
	}
	
	for (int j=0;j<6;j++) {
		mean_global_c[j] = mean_global_c[j]/streamLineList->GetNumberOfCells();
	}
	
	// En realidad min y max son el min y el max de las medias de todas las fibras calculadas
	char st[200]; 
	if (save_stats_file) Dfile.open(m_filename, ios::app);
	char name[6][4] = {"Dxx","Dxy","Dxz","Dyy","Dyz","Dzz"};
	for (int j=0;j<6;j++) {
	  sprintf(st,"Mean %f\nMax %f\nMin %f \nTot %f",mean_global_c[j],max_c[j],min_c[j],tot_c[j]);
	  m_buff[13+j]->text(st);
	  if (save_stats_file) {
		  Dfile<< name[j] << "\n"<< st << "\n------------" << std::endl;
	  }
	}
	if (save_stats_file) Dfile.close();
}


void TensorConsole::ComputeFibersFAandProject(vtkPolyData* streamLineList) {
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float mean_FA_tot, mean_FA, max_FA = 0, min_FA = 10, mean_global_FA = 0, total_FA, FA;
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	TensorImageType::RegionType region3D=aux->GetLargestPossibleRegion();
	
	FloatWriter2DType::Pointer writer = FloatWriter2DType::New();
	FloatImageType2D::Pointer ImageProjection = FloatImageType2D::New();
	FloatImageType2D::IndexType index2D;
	FloatImageType2D::SizeType size;
	size[0] = region3D.GetSize()[1];
	size[1] = region3D.GetSize()[2];
	ImageProjection->SetRegions(size);
	ImageProjection->Allocate();
	ImageProjection->FillBuffer(0);
	
	/// Para cada fibra:
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		mean_FA_tot = 0;
		//std::cout << "Numero de puntos en la siguente fibra: " << points->GetNumberOfPoints() << std::endl;
		/// Para cada punto en la fibra:
		for (int k=points->GetNumberOfPoints()-1;k>=0;k--) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			mean_FA = (float)aux->GetPixel(index).GetFractionalAnisotropy();
			index2D[0] = index[1];
			index2D[1] = index[2];
			FA = ImageProjection->GetPixel(index2D);
			ImageProjection->SetPixel(index2D,mean_FA+FA);
			mean_FA_tot=mean_FA_tot+mean_FA;
		}
		mean_FA = mean_FA_tot/(points->GetNumberOfPoints());
		if (mean_FA > max_FA) max_FA = mean_FA;
		if (mean_FA < min_FA) min_FA = mean_FA;
		mean_global_FA += mean_FA;
	}
    writer->SetInput(ImageProjection);
	writer->SetFileName("fibrasProjected.mhd");
	writer->Update();
	
	total_FA = mean_global_FA;
	mean_global_FA = mean_global_FA/streamLineList->GetNumberOfCells();
	std::cout << "-- Num of fibers computed: " << streamLineList->GetNumberOfCells()  << std::endl;
	std::cout << "Mean FA: "<< mean_global_FA << std::endl;
	std::cout << "Max FA: " << max_FA<< std::endl;
	std::cout << "Min FA: " << min_FA << std::endl;
	std::cout << "Tot FA: " << total_FA << std::endl;
	
}

void TensorConsole::SetScalarsSize(vtkPolyData* streamLineList) {
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	std::vector<float> vec_sizes;
	float dist = 0, max_size = 0, min_size = 1000000;
	double punto[3],punto_ant[3];
	
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		dist = 0;
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			if (k > 0 && k<points->GetNumberOfPoints()) {
				points->GetPoint(k,punto);
				points->GetPoint(k-1,punto_ant);
				dist += sqrt((punto_ant[0] - punto[0])*(punto_ant[0] - punto[0]) + 
			  		         (punto_ant[1] - punto[1])*(punto_ant[1] - punto[1]) +
							 (punto_ant[2] - punto[2])*(punto_ant[2] - punto[2])); 
				
			}
			
		}
		vec_sizes.push_back(dist);
		if (dist > max_size) max_size = dist;
		if (dist < min_size) min_size = dist;
	}
	
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		//std::cout << i << " Size " << vec_sizes[i] << std::endl;
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			newScalars->InsertNextValue((vec_sizes[i]));
		}
	}
	//std::cout << "NumberOfTuples "<< streamLineList->GetPointData()->GetNumberOfTuples() << std::endl;
	//std::cout << "NumberOfScalars "<< newScalars->GetNumberOfTuples() << std::endl;
	//std::cout << "min_size "<< min_size << " max_size "<< max_size << std::endl;
	if (newScalars->GetNumberOfTuples() > 0 && newScalars->GetNumberOfTuples() == streamLineList->GetPointData()->GetNumberOfTuples()) {
		//std::cout << "Setting NewScalars" << std::endl;
		streamLineList->GetPointData()->SetScalars(newScalars);
		int dataId = m_modeldataBrowser->value()-1;
		if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
		(*m_VectorModelData)[dataId].actor->GetMapper()->SetScalarRange(min_size,max_size);
		min_scalar_range->range((int)min_size,(int)max_size);
		max_scalar_range->range((int)min_size,(int)max_size);
		min_scalar_range->step((int)((int)max_size-(int)min_size)/100.0);
		max_scalar_range->step((int)((int)max_size-(int)min_size)/100.0);
		min_scalar_range->value((int)min_size);
		max_scalar_range->value((int)max_size);
		min_scalar_range->redraw();
		max_scalar_range->redraw();
	}
	
}

void TensorConsole::SetScalarsDistance(vtkPolyData* streamLineList) {
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float dist = 0, max_dist = 0;;
	double punto[3],punto_ant[3];
	
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		dist = 0;
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			if (k > 0 && k<points->GetNumberOfPoints()) {
				points->GetPoint(k,punto);
				points->GetPoint(k-1,punto_ant);
				dist += sqrt((punto_ant[0] - punto[0])*(punto_ant[0] - punto[0]) + 
			  		         (punto_ant[1] - punto[1])*(punto_ant[1] - punto[1]) +
							 (punto_ant[2] - punto[2])*(punto_ant[2] - punto[2])); 
			}
			
			if (dist > max_dist) max_dist = dist;
			newScalars->InsertNextValue(dist);	
		}
		
		
	}
	
	//for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
	//	polyline = (vtkPolyLine*)streamLineList->GetCell(i);
	//	points = polyline->GetPoints();
	//	for (int k=0;k<points->GetNumberOfPoints();k++) {
	//		newScalars->InsertNextValue((vec_sizes[i]));
	//	}
	//}
	
	std::cout << "NumberOfTuples "<< streamLineList->GetPointData()->GetNumberOfTuples() << std::endl;
	std::cout << "NumberOfScalars "<< newScalars->GetNumberOfTuples() << std::endl;
	std::cout << " max_dist "<< max_dist << std::endl;
	if (newScalars->GetNumberOfTuples() > 0 && newScalars->GetNumberOfTuples() == streamLineList->GetPointData()->GetNumberOfTuples()) {
		std::cout << "Setting Scalars Disctanc. Max: " << max_dist << std::endl;
		streamLineList->GetPointData()->SetScalars(newScalars);
		int dataId = m_modeldataBrowser->value()-1;
		if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
		(*m_VectorModelData)[dataId].actor->GetMapper()->SetScalarRange(0,max_dist);
		min_scalar_range->range(0,(int)max_dist);
		max_scalar_range->range(0,(int)max_dist);
		min_scalar_range->step((int)((int)max_dist)/100.0);
		max_scalar_range->step((int)((int)max_dist)/100.0);
		min_scalar_range->value(0);
		max_scalar_range->value((int)max_dist);
		min_scalar_range->redraw();
		max_scalar_range->redraw();
	}
	Fl::check();
}


void TensorConsole::SetScalarsFA(vtkPolyData* streamLineList) {
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			newScalars->InsertNextValue((float)aux->GetPixel(index).GetFractionalAnisotropy());
		}
		
	}
	//std::cout << "NumberOfTuples "<< streamLineList->GetPointData()->GetNumberOfTuples() << std::endl;
	//std::cout << "NumberOfScalars "<< newScalars->GetNumberOfTuples() << std::endl;
	
	if (newScalars->GetNumberOfTuples() > 0) {
		streamLineList->GetPointData()->SetScalars(newScalars);
		int dataId = m_modeldataBrowser->value()-1;
		if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
		(*m_VectorModelData)[dataId].actor->GetMapper()->SetScalarRange(0.2,0.8);
		min_scalar_range->range(0,1);
		max_scalar_range->range(0,1);
		max_scalar_range->step(0.01);
		min_scalar_range->step(0.01);
		min_scalar_range->value(0.2);
		max_scalar_range->value(0.8);
		min_scalar_range->redraw();
		max_scalar_range->redraw();
	}
	Fl::check();
}

void TensorConsole::SetScalarsMD(vtkPolyData* streamLineList) {
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			newScalars->InsertNextValue((float)aux->GetPixel(index).GetMeanDiffusivity());
		}
		
	}
	//std::cout << "NumberOfTuples "<< streamLineList->GetPointData()->GetNumberOfTuples() << std::endl;
	//std::cout << "NumberOfScalars "<< newScalars->GetNumberOfTuples() << std::endl;
	
	if (newScalars->GetNumberOfTuples() > 0) {
		streamLineList->GetPointData()->SetScalars(newScalars);
		int dataId = m_modeldataBrowser->value()-1;
		if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
		(*m_VectorModelData)[dataId].actor->GetMapper()->SetScalarRange(0.06,0.1);
		min_scalar_range->range(0.06,0.1);
		max_scalar_range->range(0.06,0.1);
		max_scalar_range->step(0.001);
		min_scalar_range->step(0.001);
		min_scalar_range->value(0.06);
		max_scalar_range->value(0.1);
		min_scalar_range->redraw();
		max_scalar_range->redraw();
	}
	Fl::check();
}

void TensorConsole::SetScalarsMajorEig(vtkPolyData* streamLineList) {
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float eig_max = 0;
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	TensorPixelType::EigenValuesArrayType eig;
	
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			aux->GetPixel(index).ComputeEigenValues(eig);
			if (eig[0] > eig_max) eig_max = eig[0];
			newScalars->InsertNextValue((float)eig[0]);
		}
		
	}
	
	if (newScalars->GetNumberOfTuples() > 0) {
		//std::cout << "Setting Scalars MajorEig. Max: " << eig_max << std::endl;
		streamLineList->GetPointData()->SetScalars(newScalars);
		int dataId = m_modeldataBrowser->value()-1;
		if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
		(*m_VectorModelData)[dataId].actor->GetMapper()->SetScalarRange(0,eig_max);
		min_scalar_range->range(0,eig_max);
		max_scalar_range->range(0,eig_max);
		max_scalar_range->step(0.01);
		min_scalar_range->step(0.01);
		min_scalar_range->value(0);
		max_scalar_range->value(eig_max);
		min_scalar_range->redraw();
		max_scalar_range->redraw();
	}
	Fl::check();
}

void TensorConsole::SetScalarsCurvature(vtkPolyData* streamLineList) {
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	float mod_vn,mod_vn_ant,c_aux,curvature;
	itk::Vector<float,3> vn,vn_ant;
	double punto[3],punto_ant[3],punto_post[3];
	
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			if (k > 0 && k<points->GetNumberOfPoints()) {
				points->GetPoint(k,punto);
				points->GetPoint(k-1,punto_ant);
				
				vn_ant[0] = punto[0] - punto_ant[0];
				vn_ant[1] = punto[1] - punto_ant[1];
				vn_ant[2] = punto[2] - punto_ant[2];
			}
			if (k> 0 && k<points->GetNumberOfPoints()-1) {
				points->GetPoint(k+1,punto_post);
			    vn[0] = punto_post[0] - punto[0];
				vn[1] = punto_post[1] - punto[1];
				vn[2] = punto_post[2] - punto[2];
			    c_aux=fabs(vn*vn_ant);
			    mod_vn=sqrt(vn*vn);
			    mod_vn_ant=sqrt(vn_ant*vn_ant);				
			    c_aux=c_aux/(mod_vn*mod_vn_ant);
			    if(c_aux>1) c_aux = 1;
			    curvature=fabs(acos( c_aux ));
				newScalars->InsertNextValue(curvature);
			    //std::cout << "Punto: " << punto[0] << " " << punto[1] << " " << punto[2] << std::endl;
			    //std::cout << "Angle: " << curvature << std::endl;
			}
		}
	}

	if (newScalars->GetNumberOfTuples() > 0) {
		streamLineList->GetPointData()->SetScalars(newScalars);
		int dataId = m_modeldataBrowser->value()-1;
		if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
		(*m_VectorModelData)[dataId].actor->GetMapper()->SetScalarRange(0,1);
	}
}

void TensorConsole::SetFibersColor(int value) {
	int dataId = m_modeldataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorModelData->size() ){return;}
	vtkPolyData* streamLineList = (*m_VectorModelData)[dataId].data;
	vtkActor* Actor = (*m_VectorModelData)[dataId].actor;
	
	
	switch(value) {
		case 0: // FA
		  SetScalarsFA(streamLineList);
		  Actor->GetMapper()->ScalarVisibilityOn();
		  break;
		case 1: // Size
		  SetScalarsSize(streamLineList);
		  Actor->GetMapper()->ScalarVisibilityOn();
		  break;
		case 2: // RGB
		  //Actor->GetMapper()->ScalarVisibilityOff();
		  //Actor->GetProperty()->SetColor(m_color_r->value(),m_color_g->value(),m_color_b->value());
		  Actor->GetMapper()->ScalarVisibilityOn();
		  break;
		case 3: // MD
		  SetScalarsMD(streamLineList);
		  Actor->GetMapper()->ScalarVisibilityOn();
	      break;
		case 4: // Distancia 
		  SetScalarsDistance(streamLineList);
		  Actor->GetMapper()->ScalarVisibilityOn();
		  break;
		case 5: // Major Eig 
		  SetScalarsMajorEig(streamLineList);
		  Actor->GetMapper()->ScalarVisibilityOn();
		  break;
		case 6:
		  Actor->GetMapper()->ScalarVisibilityOff();
		  Actor->GetProperty()->SetColor(1.0,0.8,0.2); // Skin color
		  break;
		case 7:
		  Actor->GetMapper()->ScalarVisibilityOff();
		  Actor->GetProperty()->SetColor(0.8,0.8,0.8); // bone color
		  break;
		case 8:
		  Actor->GetMapper()->ScalarVisibilityOff();
		  Actor->GetProperty()->SetColor(0.3,0.3,0.3); // gray color
		  break;
		case 9:
		  Actor->GetMapper()->ScalarVisibilityOff();
		  Actor->GetProperty()->SetColor(0.9,0.1,0.1); // red color
		  break;
	    case 10:
		  Actor->GetMapper()->ScalarVisibilityOff();
		  Actor->GetProperty()->SetColor(0.1,0.8,0.1); // green color
		  break;
		case 11:
		  Actor->GetMapper()->ScalarVisibilityOff();
		  Actor->GetProperty()->SetColor(0.1,0.5,0.7); // blue color
		  break;
		case 12:
		  Actor->GetMapper()->ScalarVisibilityOff();
		  Actor->GetProperty()->SetColor(0.1,0.1,0.7); // sea color
		  break;
		case 13: // curvature
		  SetScalarsCurvature(streamLineList);
		  Actor->GetMapper()->ScalarVisibilityOn();
		  break;
		default:
		  break;
	}
	
}

// Aplica Runge-Kutta a partir de las semillas 
void TensorConsole::RungeKuttaTractography( void ) 
{
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	std::cout << "dataId " << dataId << std::endl;
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	m_Tractography=TractographyType::New();
	//std::cout << "Doing tractography " << std::endl;
	//clicked_points = ImageViewer[m_activeinput]->getClickedPointsStored();
	InitializePoints( );
	std::list< ClickPoint >::const_iterator point = clicked_points.begin();
	
	if(!m_TractRegionsImage){
		m_TractRegionsImage=UCharImageType::New();
		m_TractRegionsImage->SetRegions(aux->GetRequestedRegion());
		m_TractRegionsImage->SetOrigin(aux->GetOrigin());
		m_TractRegionsImage->SetSpacing(aux->GetSpacing());
		m_TractRegionsImage->Allocate();
	}
	
	//std::cout << "Creating streamline " << std::endl;
	int numberOfPoints = 0;
	NodeType node;
	vtkPolyData* streamlineList=vtkPolyData::New();
	streamlineList->Allocate();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	TractographyType::InputIndexType seed;
	TractographyType::StreamlineType streamlineITK;
	
	m_Tractography->SetStepLength(m_step_length->value());
	m_Tractography->SetFaThreshold(m_fa_threshold->value());
	m_Tractography->SetCurvatureThreshold(m_curvature_threshold->value());
	m_Tractography->SetInput(aux);
	
	//std::cout << "Haciendo tractografia del volumen con dataId: " << dataId << " cuyo Origen de coordenadas es " << aux->GetOrigin() << std::endl;
	
	while( point != clicked_points.end() ) {
		seed[0]=round(point->x);
		seed[1]=round(point->y);
		seed[2]=round(point->z);
		
		vtkIdList* streamline=vtkIdList::New();
		streamline->Allocate(0);
		
		m_Tractography->RungeKuttaIntegration(seed, streamlineITK);
		TensorImageType::SpacingType spacing=aux->GetSpacing();
		TensorImageType::PointType origen=aux->GetOrigin();
		if (streamlineITK.size() > 2) { 
			for ( unsigned j = 0; j < streamlineITK.size(); j++ ) {
				points->InsertNextPoint(streamlineITK[j][0],streamlineITK[j][1],streamlineITK[j][2]);
				TensorImageType::IndexType index;
				index[0] = round((streamlineITK[j][0]-origen[0])/spacing[0]);
				index[1] = round((streamlineITK[j][1]-origen[1])/spacing[1]);
				index[2] = round((streamlineITK[j][2]-origen[2])/spacing[2]);
				TensorImageType::SizeType Size = aux->GetLargestPossibleRegion().GetSize();
				if (colortype[0]->value()) {
					// Si quiero color por FA:
					if (index[0] >= 0 && index[0] < (int)Size[0] && 
						index[1] >= 0 && index[1] < (int)Size[1] &&
						index[2] >= 0 && index[2] < (int)Size[2]) {
						newScalars->InsertNextValue((float)aux->GetPixel(index).GetFractionalAnisotropy());
					}
				}
				if (colortype[1]->value()) {
					//Si quiero color por tamaño:
					//std::cout << "fibra " << j << " tamaño: " << streamlineITK.size()/100 << std::endl;
					newScalars->InsertNextValue(float(streamlineITK.size())/100);
				}
				if (colortype[2]->value()) {
					//Si quiero color arbitrario:
					//std::cout << "fibra " << j << " valor: " << (float)m_imageOverlay->GetPixel(seed) << std::endl;
					newScalars->InsertNextValue((float)m_imageOverlay->GetPixel(seed)/30);
				}
				
				if (colortype[3]->value()) {
					//Si quiero color por curvatura:
					std::cout << "fibra " << j << " curvatura: " << streamlineITK.size() << std::endl;
					newScalars->InsertNextValue(float(streamlineITK.size()));
				}
				if (colortype[4]->value()) {
					//Si quiero color por distancia a la semilla:
					//std::cout << "fibra " << j << " distancia: " << float(j)/float(streamlineITK.size()) << std::endl;
					newScalars->InsertNextValue(float(j)/float(streamlineITK.size()));
				}
				if (colortype[5]->value()) {
					//Si quiero color por mayor autovalor:
					if (index[0] >= 0 && index[0] < (int)Size[0] && 
						index[1] >= 0 && index[1] < (int)Size[1] &&
						index[2] >= 0 && index[2] < (int)Size[2]) {
						TensorPixelType::EigenValuesArrayType eig;
						aux->GetPixel(index).ComputeEigenValues(eig);
						//std::cout << "punto " << j << " major eigenvalue: " << (float)eig[0]*100 << std::endl;
						newScalars->InsertNextValue((float)eig[0]*100);
					}
				}
				if (index[0] >= 0 && index[0] < (int)Size[0] && 
					index[1] >= 0 && index[1] < (int)Size[1] &&
					index[2] >= 0 && index[2] < (int)Size[2]) {
					m_TractRegionsImage->SetPixel(index, m_VectorModelData->size()+1);
				}
				streamline->InsertNextId(numberOfPoints);
				numberOfPoints++;
			}
			
			streamlineList->InsertNextCell(VTK_POLY_LINE, streamline);
		}
		streamline->Delete();
		streamlineITK.clear();
		
		++point;		
	} // end while
	//std::cout << "Num de puntos totales : " << numberOfPoints << std::endl;

	///////// 
	streamlineList->SetPoints(points);
	//ComputeFibersSize(streamlineList);
	//ComputeFibersAngle(streamlineList);
	//ComputeFibersFA(streamlineList,false);
	
	if (colortype[0]->value() || colortype[1]->value() || colortype[2]->value() || colortype[3]->value() || colortype[4]->value() || colortype[5]->value() ) {
		streamlineList->GetPointData()->SetScalars(newScalars);
	}
	
	char nombre[200];
	sprintf(nombre,"Streamline %d",(int)m_VectorModelData->size()); 
	
	vtkActor *Actor = vtkActor::New();
	DataModelElementType model_data;
	model_data.Id     = m_VectorModelData->size();
	model_data.nombre = nombre;
	model_data.data   = streamlineList;
	model_data.actor  = Actor;
	
	m_VectorModelData->push_back(model_data);
	
	GetColorFromSeed((unsigned int)m_imageOverlay->GetPixel(seed));
	double color[3] = {m_color_r->value(),m_color_g->value(),m_color_b->value()};
	//std::cout << (unsigned int)m_imageOverlay->GetPixel(seed) << " color " << color[0] << " "<< color[1] << " " << color[2] << std::endl;	
	
	ImageViewer3D->renderTracts(streamlineList, Actor, color, m_radius->value());

	m_puntosTract = vtkPoints::New();
	m_puntosTract->DeepCopy(points);
	
	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();
}


void TensorConsole::BruteForceTractography() 
{
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	
	m_Tractography=TractographyType::New();
	std::cout << "Doing tractography all" << std::endl;
	m_Tractography->SetNumberOfThreads(1);
	m_Tractography->SetStepLength(m_step_length->value());
	m_Tractography->SetFaThreshold(m_fa_threshold->value());
	m_Tractography->SetCurvatureThreshold(m_curvature_threshold->value());
	m_Tractography->SetInput((*m_VectorTensorData)[dataId].image);
	
	try{
		m_Tractography->Update();
		m_HasBeenFiberTracked=true;
	}
	catch ( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
}

void TensorConsole::SelectFiberTracts(int label) 
{
	if(!m_Tractography){
		return;
	}
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	
	if(!m_TractRegionsImage){
		m_TractRegionsImage=UCharImageType::New();
		m_TractRegionsImage->SetRegions(aux->GetRequestedRegion());
		m_TractRegionsImage->SetOrigin(aux->GetOrigin());
		m_TractRegionsImage->SetSpacing(aux->GetSpacing());
		m_TractRegionsImage->Allocate();
	}
	
	
	//std::list<ClickPoint> clicked_points;
	//clicked_points = ImageViewer[m_activeinput]->getClickedPointsStored();
	if (label < 0) { 
		InitializePoints(  );
	} else {
 	    InitializePoints( label );
	}
	std::list< ClickPoint >::const_iterator point = clicked_points.begin();
	
	int numberOfPoints = 0;
	NodeType node;
	vtkPolyData* streamlineList=vtkPolyData::New();
	streamlineList->Allocate();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	TractographyType::InputIndexType seed;
	TractographyType::StreamlineType streamlineITK;
	
	
	VectorImageType::Pointer fiberIndexImage=m_Tractography->GetOutput();
	typedef itk::ImageRegionConstIterator<VectorImageType> VectorConstIteratorType;
	
	TractographyType::StreamlineListType	streamlineImage;
	m_Tractography->GetStreamlineList(streamlineImage);		
	
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	TensorImageType::PointType pointAct;
	TensorImageType::PointType pointAnt;
	
	std::vector<float> vec_sizes;
	float mean_size = 0;
	char filename[25];
	sprintf(filename, "longitudFibras%d.txt", (label-1));
	std::cout<<filename<<std::endl; 
	ofstream strm(filename);
	
	while( point != clicked_points.end() ) {
		seed[0]=round(point->x);
		seed[1]=round(point->y);
		seed[2]=round(point->z);
		itk::VariableLengthVector<unsigned long> streamlineIndex;
		streamlineIndex=fiberIndexImage->GetPixel(seed) ;
		unsigned i=0;
		float dist=0;
		
		if(streamlineIndex.GetSize()!=0){
			while(streamlineIndex[i]>0 && i < streamlineIndex.GetSize()){
				dist=0;
				vtkIdList* streamline=vtkIdList::New();
				streamline->Allocate(0);
				
				if( ( streamlineIndex[i]-1) < 0 || (streamlineIndex[i]-1) > streamlineImage.size() ){
					std::cout<<i<<" "<<streamlineIndex[i]<<std::endl;
				}else{
					streamlineITK=streamlineImage[streamlineIndex[i]-1];
					//std::cout<<"numero puntos streamline "<< streamlineITK.size() <<std::endl;
					
					for ( unsigned j = 0; j < streamlineITK.size(); j++ ) {
						points->InsertNextPoint(streamlineITK[j][0],streamlineITK[j][1],streamlineITK[j][2]);
						pointAct[0]=streamlineITK[j][0];
						pointAct[1]=streamlineITK[j][1];
						pointAct[2]=streamlineITK[j][2];
						
						if(j>0){
							dist += sqrt((pointAnt[0] - pointAct[0])*(pointAnt[0] - pointAct[0]) + 
										 (pointAnt[1] - pointAct[1])*(pointAnt[1] - pointAct[1]) +
										 (pointAnt[2] - pointAct[2])*(pointAnt[2] - pointAct[2])); 
						}
						TensorImageType::IndexType index;
						index[0] = round((streamlineITK[j][0]-origen[0])/spacing[0]);
						index[1] = round((streamlineITK[j][1]-origen[1])/spacing[1]);
						index[2] = round((streamlineITK[j][2]-origen[2])/spacing[2]);
						TensorImageType::SizeType Size = aux->GetLargestPossibleRegion().GetSize();
						if (index[0] >= 0 && index[0] < (int)Size[0] && 
							index[1] >= 0 && index[1] < (int)Size[1] &&
							index[2] >= 0 && index[2] < (int)Size[2]) {
							m_TractRegionsImage->SetPixel(index, label);
						}
						
						if (colortype[0]->value()) {
							// Si quiero color por FA:
							if (index[0] >= 0 && index[0] < (int)Size[0] && 
								index[1] >= 0 && index[1] < (int)Size[1] &&
								index[2] >= 0 && index[2] < (int)Size[2]) {
								newScalars->InsertNextValue((float)aux->GetPixel(index).GetFractionalAnisotropy());
								
							}
						}
						if (colortype[1]->value()){
							//Si quiero color por tamaño:
							newScalars->InsertNextValue(float(j)/float(streamlineITK.size()));
						}
						if (colortype[2]->value()) {
							//Si quiero color arbitrario:
							//std::cout << "fibra " << j << " valor: " << (float)m_imageOverlay->GetPixel(seed) << std::endl;
							newScalars->InsertNextValue((float)m_imageOverlay->GetPixel(seed)/31);
						}
						
						streamline->InsertNextId(numberOfPoints);
						numberOfPoints++;
						pointAnt=pointAct;
					}					
					if (streamline->GetNumberOfIds() > 0) {	
						streamlineList->InsertNextCell(VTK_POLY_LINE, streamline);
					}		
					streamline->Delete();
					vec_sizes.push_back(dist);
					mean_size+=dist;
				}
				
				i++;
			}
		}
		++point;		
	} // end while
	
	//	std::cout << "Num puntos: " << numberOfPoints << std::endl;
	
	///////// 
	streamlineList->SetPoints(points);
	if (colortype[0]->value() || colortype[1]->value() || colortype[2]->value() || colortype[3]->value() || colortype[4]->value() || colortype[5]->value() ) {
		streamlineList->GetPointData()->SetScalars(newScalars);
	}
	
	mean_size=mean_size/streamlineList->GetNumberOfCells();
	std::cout<<"tamaño de fibra medio: "<<mean_size<<std::endl;
	strm<<"tracto: "<<label-1<<" tamano :"<<mean_size<<std::endl;
	
	char nombre[200];
	sprintf(nombre,"Streamline %d",(int)m_VectorModelData->size());
	
	vtkActor *Actor = vtkActor::New();
	DataModelElementType model_data;
	model_data.Id     = m_VectorModelData->size();
	model_data.nombre = nombre;
	model_data.data   = streamlineList;
	model_data.actor  = Actor;
	m_VectorModelData->push_back(model_data);
	
	/**
	 m_modeldataBrowser->add(nombre);
	 m_modeldataBrowser->select(m_modeldataBrowser->size());
	 m_modeldataBrowser->redraw();
	 */
	GetColorFromSeed(m_imageOverlay->GetPixel(seed));
	double color[3] = {m_color_r->value(),m_color_g->value(),m_color_b->value()};
	ImageViewer3D->renderTracts(streamlineList, Actor, color, m_radius->value());
	
	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();
}

void TensorConsole::SelectAutoFiberTracts(int label1) { 

	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	
	if(!m_HasBeenFiberTracked){
		try
		{
			BruteForceTractography();
			m_HasBeenFiberTracked=true;
		}
		catch (itk::ExceptionObject & e)
		{
			fl_alert( e.GetDescription() );
			return;
		}
	}
	
	std::vector<int> rois_excluded; 
	
	switch (label1) {
		case 1: 
			rois_excluded.push_back(29);
			SelectFiberTractsConnection(label1,32,rois_excluded);
			break;
		case 15:
			rois_excluded.push_back(29);
			SelectFiberTractsConnection(label1,33,rois_excluded);
			break;
		case 12: // Fasciculo uncinado izquierdo
			rois_excluded.push_back(9);
			rois_excluded.push_back(29);
			SelectFiberTractsConnection(label1,label1,rois_excluded);
		    break;
		case 26: // fasciculo uncinado derecho
			rois_excluded.push_back(23);
			rois_excluded.push_back(29);
			SelectFiberTractsConnection(label1,label1,rois_excluded);			
			break;
		case 9: // ifo izquierdo
			rois_excluded.push_back(10);
			rois_excluded.push_back(29);
			SelectFiberTractsConnection(label1,label1,rois_excluded);			
			break;
		case 23: // ifo dcho
			rois_excluded.push_back(24);
			rois_excluded.push_back(29);
			SelectFiberTractsConnection(label1,label1,rois_excluded);			
			break;
		case 10: // Slf izquierdo
			rois_excluded.push_back(9);
			rois_excluded.push_back(29);
			SelectFiberTractsConnection(label1,label1,rois_excluded);			
			break;
		case 24: // Slf derecho
			rois_excluded.push_back(23);
			rois_excluded.push_back(29);
			SelectFiberTractsConnection(label1,label1,rois_excluded);			
			break;
		case 29: // Cuerpo Calloso
			rois_excluded.push_back(1);
			rois_excluded.push_back(15);
			SelectFiberTractsConnection(label1,label1,rois_excluded);			
			break;
		default:
			SelectFiberTracts(label1);
			break;
	}
	
	
}


void TensorConsole::SelectFiberTractsConnection(int label1, int label2, std::vector<int> rois_excluded) 
{
	if(!m_Tractography){
		return;
	}
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	
	InitializePoints( );
	std::list< ClickPoint >::const_iterator point = clicked_points.begin();
	
	int numberOfPoints = 0;
	NodeType node;
	vtkPolyData* streamlineList=vtkPolyData::New();
	streamlineList->Allocate();
	//vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	vtkFloatArray *newScalars = vtkFloatArray::New();
	newScalars->Allocate(0);
	TractographyType::InputIndexType seed;
	TractographyType::StreamlineType streamlineITK;
	
	VectorImageType::Pointer fiberIndexImage=m_Tractography->GetOutput();
	TractographyType::StreamlineListType	streamlineImage;
	m_Tractography->GetStreamlineList(streamlineImage);		
	
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	TensorImageType::PointType pointAct;
	TensorImageType::PointType pointAnt;
	
	std::vector<unsigned long> vectorA;
	std::vector<unsigned long> vectorB;
	std::vector<unsigned long> vectorExcluded;
	
	//std::cout << "rois_excluded: size " <<  rois_excluded.size() << " value[0] " << rois_excluded[0] << std::endl;

	
	while( point != clicked_points.end() ) {
		seed[0]=round(point->x);
		seed[1]=round(point->y);
		seed[2]=round(point->z);
		itk::VariableLengthVector<unsigned long> streamlineIndex;
		
		if (m_imageOverlay->GetPixel(seed) == label1) {
			streamlineIndex=fiberIndexImage->GetPixel(seed);
			if (vectorA.size() == 0) {
				vectorA.push_back(streamlineIndex[0]);
			}
			for (unsigned int i=0;i<streamlineIndex.GetSize();i++) {
				bool flag = true;
				for (unsigned int j=0;j<vectorA.size();j++) {
					if (vectorA[j] == streamlineIndex[i]) {
						flag = false;
						break;
					}
			    }
				if (flag) vectorA.push_back(streamlineIndex[i]); 
			}
		}
		
		if (m_imageOverlay->GetPixel(seed) == label2) {
			streamlineIndex=fiberIndexImage->GetPixel(seed);
			if (vectorB.size() == 0) {
				vectorB.push_back(streamlineIndex[0]);
			}
			for (unsigned int i=0;i<streamlineIndex.GetSize();i++) {
				bool flag = true;
				for (unsigned int j=0;j<vectorB.size();j++) {
					if (vectorB[j] == streamlineIndex[i]) {
						flag = false;
						break;
					}
			    }
				if (flag) vectorB.push_back(streamlineIndex[i]); 
			}		
		}
		
		if (rois_excluded.size() > 0) {
		  for (unsigned int k=0;k<rois_excluded.size();k++) {
		    if (m_imageOverlay->GetPixel(seed) == rois_excluded[k] ) {
			  streamlineIndex=fiberIndexImage->GetPixel(seed);
			  if (vectorExcluded.size() == 0) {
				vectorExcluded.push_back(streamlineIndex[0]);
			  }
			  for (unsigned int i=0;i<streamlineIndex.GetSize();i++) {
				bool flag = true;
				for (unsigned int j=0;j<vectorExcluded.size();j++) {
				  if (vectorExcluded[j] == streamlineIndex[i]) {
					flag = false;
					break;
				  }
			    } // end for j
				if (flag) vectorExcluded.push_back(streamlineIndex[i]); 
			  } // end for i
		    } // end if m_imageOverlay->GetPixel(seed) == rois_excluded[k]
		  }  // end for k
		} //end if rois_excluded.size() > 0
		++point;
	}
	
	//std::cout << "vectorA: size " <<  vectorA.size() << std::endl;
	//	for (unsigned int i=0;i< vectorA.size();i++) std::cout << vectorA[i] << " "; 	
	//	std::cout << std::endl;
	//	
	//	
	//	std::cout << "vectorB: size " <<  vectorB.size() <<  std::endl;
	//	for (unsigned int i=0;i< vectorB.size();i++) std::cout << vectorB[i] << " "; 	
	std::cout << std::endl;
	
	//Juntemos ahora las que coinciden
	unsigned int s=0;
	std::vector<unsigned long> vectorOut;
	for (unsigned int i=0;i<vectorA.size();i++) {
		for (unsigned int j=0;j<vectorB.size();j++) {
			if (vectorB[j] == vectorA[i]) {
				vectorOut.push_back(vectorA[i]);
				s++;
				break;
			}
		}	
	}
	std::cout << "vectorOut: size " <<  vectorOut.size() << " s " << s << std::endl;
	//	for (unsigned int i=0;i< vectorOut.size();i++) std::cout << vectorOut[i] << " "; 	
	//	std::cout << std::endl;
	
	//Quitamos ahora las que coinciden con las regiones que queremos excluir
	s=0;
	//std::vector<unsigned long> vectorAux = vectorOut;
	std::vector<unsigned long>::iterator It; 
	for (unsigned int i=0;i<vectorExcluded.size();i++) {
	  It = vectorOut.begin();
	  for (unsigned int j=0;j<vectorOut.size();j++) {
		//while ( It != vectorOut.end() ) {
			if (vectorExcluded[i] == vectorOut[j]) {
				vectorOut.erase(It);
				s++;
				break;
			}
			It++;
		} 	
	}
	std::cout << "vectorOut: size " <<  vectorOut.size() << " s " << s << std::endl;
	
	float mean_size = 0, min_size = 10000000, max_size = 0;
	std::vector<float> vec_sizes;
	
	//Calculo el tamaño medio de todas las fibras
	for (unsigned int i=0;i<vectorOut.size();i++){
		float dist = 0;
		
		if( ( vectorOut[i]-1) < 0 || (vectorOut[i]-1) > streamlineImage.size() ){
			std::cout << i <<" "<< vectorOut[i]<<std::endl;
		}else{
			streamlineITK=streamlineImage[vectorOut[i]-1];
			
			for ( unsigned j = 0; j < streamlineITK.size(); j++ ) {
				pointAct[0]=streamlineITK[j][0];
				pointAct[1]=streamlineITK[j][1];
				pointAct[2]=streamlineITK[j][2];
				
				if(j>0){
					dist += sqrt((pointAnt[0] - pointAct[0])*(pointAnt[0] - pointAct[0]) + 
				  		         (pointAnt[1] - pointAct[1])*(pointAnt[1] - pointAct[1]) +
								 (pointAnt[2] - pointAct[2])*(pointAnt[2] - pointAct[2])); 
				}	
				pointAnt=pointAct;
			}					
			vec_sizes.push_back(dist);
			
			mean_size+=dist;
			if (dist > max_size) max_size = dist;
			if (dist < min_size) min_size = dist;
		}
	}
	mean_size = mean_size/(vectorOut.size()-1);	
	
	// 	std::cout << "-- Num of fibers computed: " << vectorOut.size()-1 << std::endl;
	// 	std::cout << "Mean size of all fibers: " << mean_size << std::endl;
	// 	std::cout << "Max size of all fibers: " << max_size << std::endl;
	// 	std::cout << "Min size of all fibers: " << min_size << std::endl;
	
	
	// Ahora metemos las fibras que quiero en vtk para pintarlas
	for (unsigned int i=1;i<vectorOut.size();i++){
		if (vec_sizes[i-1]<1.4*mean_size){
			vtkIdList* streamline=vtkIdList::New();
			streamline->Allocate(0);
			
			if( ( vectorOut[i]-1) < 0 || (vectorOut[i]-1) > streamlineImage.size() ){
				std::cout << i <<" "<< vectorOut[i]<<std::endl;
			}else{
				streamlineITK=streamlineImage[vectorOut[i]-1];
				//std::cout<<"numero puntos streamline "<< streamlineITK.size() <<std::endl;
				
				for ( unsigned j = 0; j < streamlineITK.size(); j++ ) {
					points->InsertNextPoint(streamlineITK[j][0],streamlineITK[j][1],streamlineITK[j][2]);
					if (colortype[0]->value()) {
						// Si quiero color por FA:
						TensorImageType::IndexType index;
						index[0] = round((streamlineITK[j][0]-origen[0])/spacing[0]);
						index[1] = round((streamlineITK[j][1]-origen[1])/spacing[1]);
						index[2] = round((streamlineITK[j][2]-origen[2])/spacing[2]);
						TensorImageType::SizeType Size = aux->GetLargestPossibleRegion().GetSize();
						if (index[0] >= 0 && index[0] < (int)Size[0] && 
							index[1] >= 0 && index[1] < (int)Size[1] &&
							index[2] >= 0 && index[2] < (int)Size[2]) {
							newScalars->InsertNextValue((float)aux->GetPixel(index).GetFractionalAnisotropy());
						}
					}
					if (colortype[1]->value()){
						//Si quiero color por tamaño:
						newScalars->InsertNextValue(float(j)/float(streamlineITK.size()));
					}
					if (colortype[2]->value()) {
						//Si quiero color arbitrario:
						//newScalars->InsertNextValue((float)m_imageOverlay->GetPixel(seed)/31);
						newScalars->InsertNextValue((float)(vectorOut[i]-1)/31);
					}
					streamline->InsertNextId(numberOfPoints);
					numberOfPoints++;
				}
				if (streamline->GetNumberOfIds() > 0) {	
					streamlineList->InsertNextCell(VTK_POLY_LINE, streamline);
				}
				streamline->Delete();
			}
		}//else{std::cout<<"Tamano de fibra "<<i<<" demasiado larga: "<<vec_sizes[i-1]<<std::endl;}
	}
	
	//std::cout << "Num puntos: " << numberOfPoints << std::endl;
	
	
	streamlineList->SetPoints(points);
	///////// Para mostrar estadisticas por pantalla:
	//ComputeFibersSize(streamlineList);
	//ComputeFibersAngle(streamlineList);
	//ComputeFibersFA(streamlineList, false);
	
	
	if (colortype[0]->value() || colortype[1]->value() || colortype[2]->value() || colortype[3]->value() || colortype[4]->value() || colortype[5]->value() ) {
		streamlineList->GetPointData()->SetScalars(newScalars);
	}
	
	char nombre[200];
	sprintf(nombre,"Streamline %d",(int)m_VectorModelData->size());
	
	vtkActor *Actor = vtkActor::New();
	DataModelElementType model_data;
	model_data.Id     = m_VectorModelData->size();
	model_data.nombre = nombre;
	model_data.data   = streamlineList;
	model_data.actor  = Actor;
	m_VectorModelData->push_back(model_data);
	
	/**
	 m_modeldataBrowser->add(nombre);
	 m_modeldataBrowser->select(m_modeldataBrowser->size());
	 m_modeldataBrowser->redraw();
	 */
	
	double color[3] = {m_color_r->value(),m_color_g->value(),m_color_b->value()};
	ImageViewer3D->renderTracts(streamlineList, Actor, color, m_radius->value());
	
	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();
	
}


void TensorConsole::ColorOrientation(void)
{
	
	int dataId = m_tensordataBrowser->value()-1;
	if(dataId<0){return;}
	TensorImageType::Pointer tensor = (*m_VectorTensorData)[dataId].image;
	
	ImageRGBType::Pointer image = ImageRGBType::New();
	image->CopyInformation( tensor );
	image->SetRegions( tensor->GetLargestPossibleRegion() );
	image->SetMetaDataDictionary( tensor->GetMetaDataDictionary() );
	
	try{
		image->Allocate();
	}
	catch ( itk::ExceptionObject & e ){
		fl_alert( e.GetDescription() );
		return;
	}
	
	
	
	RGBIteratorType imageIt(image,image->GetRequestedRegion());
	ConstTensorIteratorType tensorIt(tensor,tensor->GetRequestedRegion());
	for (tensorIt.GoToBegin(),imageIt.GoToBegin();!tensorIt.IsAtEnd();++tensorIt, ++imageIt) {
		if (tensorIt.GetIndex()[0] == 10 && tensorIt.GetIndex()[1] == 10 && tensorIt.GetIndex()[2] == 10) {
			std::cout << tensorIt.Get()[0] << " " << tensorIt.Get()[1] << " " << tensorIt.Get()[2] << "\n";
			std::cout << tensorIt.Get()[1] << " " << tensorIt.Get()[3] << " " << tensorIt.Get()[4] << "\n";
			std::cout << tensorIt.Get()[2] << " " << tensorIt.Get()[4] << " " << tensorIt.Get()[5] << "\n";
			imageIt.Set(tensorIt.Get().GetRGBColorCode());
			TensorPixelType::EigenValuesArrayType eigval;
			
			TensorPixelType::EigenVectorsMatrixType eigvec;
			tensorIt.Get().ComputeEigenSystem( eigval, eigvec );
			
			std::cout << "eigval 1: " << eigval[0] << "\n";
			std::cout << "eigvec 1: " << eigvec[0][0] << " " << eigvec[1][0] << " " << eigvec[2][0] << "\n";
			std::cout << "eigval 2: " << eigval[1] << "\n";
			std::cout << "eigvec 2: " << eigvec[0][1] << " " << eigvec[1][1] << " " << eigvec[2][1] << "\n";
			std::cout << "eigval 3: " << eigval[2] << "\n";
			std::cout << "eigvec 3: " << eigvec[0][2] << " " << eigvec[1][2] << " " << eigvec[2][2] << "\n";
			std::cout << "%%%%%%%%%%%% \n ";
			TensorPixelType::RealType eigvec_[3];
			tensorIt.Get().ComputeEigenValues(eigval);
			tensorIt.Get().ComputeEigenVector( eigval[0], eigvec_);
			std::cout << "eigval 1: " << eigval[0] << "\n";
			std::cout << "eigvec 1: " << eigvec_[0] << " " << eigvec_[1] << " " << eigvec_[2] << "\n";
			tensorIt.Get().ComputeEigenVector( eigval[2], eigvec_);
			std::cout << "eigval 2: " << eigval[1] << "\n";
			std::cout << "eigvec 2: " << eigvec_[0] << " " << eigvec_[1] << " " << eigvec_[2] << "\n";
			tensorIt.Get().ComputeEigenVector( eigval[1], eigvec_);
			std::cout << "eigval 3: " << eigval[2] << "\n";
			std::cout << "eigvec 3: " << eigvec_[0] << " " << eigvec_[1] << " " << eigvec_[2] << "\n";
			
		} else {
			imageIt.Set(tensorIt.Get().GetRGBColorCode());
		}
	}
	
	ImageColorViewer->SetInputImage(image);
	std::cout << "ImageColorViewer slices: " << ImageColorViewer->numSlices() << std::endl;
	
	sliceNumberSlider[3]->range( 0.0f, ImageColorViewer->numSlices()-1 );
	sliceNumberSlider[3]->value((float)ImageColorViewer->sliceNum());
	
	ImageColorViewer->Show();
	ImageColorViewer->ViewDetails(false);
	ImageColorViewer->ViewValue(false);
	
	ImageColorViewer->update();
	Fl::check();
	ImageColorViewer->redraw();
	Fl::check();
	
}

void TensorConsole::CreateAutomaticFiberTractRuben(std::vector<unsigned int> tractIndex){
	
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	
	if(!m_HasBeenRegistered){
		try
		{
			RegisterSeedsToData();
			m_HasBeenRegistered=true;
			
		}
		catch ( itk::ExceptionObject & e )
		{
			fl_alert( e.GetDescription() );
			return;
		}
	}
	
	if(!m_HasBeenFiberTracked){
		try
		{
			BruteForceTractography();
			m_HasBeenFiberTracked=true;
		}
		catch (itk::ExceptionObject & e)
		{
			fl_alert( e.GetDescription() );
			return;
		}
	}
	
	m_imageOverlay=m_SeedRegionsImage;
	
	// grabamos las ROIs deformadas 
	typedef itk::ImageFileWriter< UCharImageType >   WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("ROIs_deform.mhd");
	writer->SetInput(m_SeedRegionsImage);
	writer->Update();
	
	// De momento solo para la conexion del tracto piramidal
	std::vector<int> fibers_excluded; 
	SelectFiberTractsConnection(1,32, fibers_excluded);
	SelectFiberTractsConnection(15,33, fibers_excluded);
	
	
	if(!m_TractRegionsImage){
		m_TractRegionsImage=UCharImageType::New();
		m_TractRegionsImage->SetRegions(m_SeedRegionsImage->GetRequestedRegion());
		m_TractRegionsImage->Allocate();
	}		
	
	
	/*	typedef itk::ImageFileWriter<UCharImageType> escribeType;
	 escribeType::Pointer escribe=escribeType::New();
	 escribe->SetFileName("regiones.vtk");
	 escribe->SetInput(m_TractRegionsImage);
	 escribe->Update();*/
	bool evaluate=false;
	bool correlate=false;
	bool angleCoherence=false;
	
	if(evaluate){
		//EvaluateTractRegion(tractIndex,0);
	}
	if(correlate){
		CorrelateTractRegion(tractIndex,0);
	}
	if(angleCoherence){
		ComputeAngleCoherenceTractRegion(tractIndex,0);
	}
	
};

void TensorConsole::CreateAutomaticFiberTract(std::vector<unsigned int> tractIndex){
	for (unsigned int i = 0;i< tractIndex.size();i++) {
  		std::cout << "tractIndex " << tractIndex[i] << std::endl;
	}
	
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	
	if(!m_HasBeenRegistered){
		try
		{
			RegisterSeedsToData();
			m_HasBeenRegistered=true;
			
		}
		catch ( itk::ExceptionObject & e )
		{
			fl_alert( e.GetDescription() );
			return;
		}
	}
	
	if(!m_HasBeenFiberTracked){
		try
		{
			BruteForceTractography();
			m_HasBeenFiberTracked=true;
		}
		catch (itk::ExceptionObject & e)
		{
			fl_alert( e.GetDescription() );
			return;
		}
	}
	
	m_imageOverlay=m_SeedRegionsImage;
	
	for(unsigned i=0; i<tractIndex.size(); ++i){
		SelectFiberTracts(tractIndex[i]+1);
	}
	
	if(!m_TractRegionsImage){
		m_TractRegionsImage=UCharImageType::New();
		m_TractRegionsImage->SetRegions(m_SeedRegionsImage->GetRequestedRegion());
		m_TractRegionsImage->Allocate();
	}		
	
	
	/*	typedef itk::ImageFileWriter<UCharImageType> escribeType;
	 escribeType::Pointer escribe=escribeType::New();
	 escribe->SetFileName("regiones.vtk");
	 escribe->SetInput(m_TractRegionsImage);
	 escribe->Update();*/
	bool evaluate=true;
	bool correlate=false;
	bool angleCoherence=false;
	
	if(evaluate){
		//EvaluateTractRegion(tractIndex,0);
	}
	if(correlate){
		CorrelateTractRegion(tractIndex,0);
	}
	if(angleCoherence){
		ComputeAngleCoherenceTractRegion(tractIndex,0);
	}
	
};

void TensorConsole::ProjectROIsFromFibers(vtkPolyData* streamLineList) 
{
	
	vtkPolyLine* polyline = vtkPolyLine::New();
	vtkPoints* points= vtkPoints::New();
	points->Allocate(0);
	double punto[3];
	TensorImageType::IndexType index;
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	
	TensorImageType::Pointer aux = (*m_VectorTensorData)[dataId].image;
	TensorImageType::SpacingType spacing=aux->GetSpacing();
	TensorImageType::PointType origen=aux->GetOrigin();
	
	unsigned int value = m_modeldataBrowser->value();
	
	for (int i=0;i<streamLineList->GetNumberOfCells();i++) {
		polyline = (vtkPolyLine*)streamLineList->GetCell(i);
		points = polyline->GetPoints();
		for (int k=0;k<points->GetNumberOfPoints();k++) {
			points->GetPoint(k,punto);
			index[0] = round((punto[0]-origen[0])/spacing[0]);
			index[1] = round((punto[1]-origen[1])/spacing[1]);
			index[2] = round((punto[2]-origen[2])/spacing[2]);
			m_imageOverlay->SetPixel(index,value);
		}
		
	}
	
	Fl::check();
	ImageViewer[0]->redraw();
	Fl::check();
}


void TensorConsole::EvaluateTractRegion(std::vector<unsigned int> tractIndex, std::vector<unsigned int> ind_measures, bool coherence){
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	
	ComputeMeasuresPointer measures=ComputeMeasuresType::New();
	measures->SetInput((*m_VectorTensorData)[dataId].image);
	
	//measures->SetInput2(m_TractRegionsImage);
	measures->SetInput2(m_imageOverlay);
	
	//itk::ImageRegionConstIterator<UCharImageType>  mit = itk::ImageRegionConstIterator<UCharImageType>( m_TractRegionsImage, m_TractRegionsImage->GetRequestedRegion() );
	itk::ImageRegionConstIterator<UCharImageType>  mit = itk::ImageRegionConstIterator<UCharImageType>( m_imageOverlay, m_imageOverlay->GetRequestedRegion() );
	itk::ImageRegionConstIterator<FloatImageType>  it;
	
	//ofstream strm1("coherenciaMedidas.txt");
	//ofstream strm2("ListaMedidas.txt");
	std::string  strm1_new;
	std::string  strm2_new;
	char s[3];
	char sf[10];
	
	strm1_new = "   ";
	strm2_new = "   ";
	for(unsigned i=0; i<tractIndex.size(); ++i){
		sprintf(s,"%02d",tractIndex[i]+1);
		strm1_new += "    " + std::string(s) + "   ";
	    strm2_new += "    " + std::string(s) + "   " ;
	}
	strm1_new += "\n";
	strm2_new += "\n";
	
	std::vector<float> avg_measure;
	for(unsigned int i=0; i<ind_measures.size(); i++){
		//if(!coherence){
			//strm1<<"Medida "<< i <<":"<<std::endl;
		//}else{
			//strm2<<"Medida "<< i<<":"<<std::endl;
		//}
		switch(ind_measures[i]){
			case 0:
				measures->SetComputeFA( );
				if (!coherence) strm1_new += "FA: ";
				if (coherence) strm2_new += "FA: ";
				break;
			case 1:
				measures->SetComputeRA( );
				if (!coherence) strm1_new +="RA: ";
				if (coherence) strm2_new +="RA: ";
				break;
			case 2:
				measures->SetComputeMD( );
				if (!coherence) strm1_new +="MD: ";
				if (coherence) strm2_new +="MD: ";
				break;
			case 3:
				measures->SetComputeEigVal( 0 );
				if (!coherence) strm1_new +="E1: ";
				if (coherence) strm2_new +="E1: ";
				break;
			case 4:
				measures->SetComputeEigVal( 1 );
				if (!coherence) strm1_new +="E2: ";
				if (coherence) strm2_new +="E2: ";
				break;
			case 5:
				measures->SetComputeEigVal( 2 );
				if (!coherence) strm1_new +="E3: ";
				if (coherence) strm2_new +="E3: ";
				break;
			case 6:
				measures->SetComputeShapeCoefficients( 0 );
				if (!coherence) strm1_new +="Cl: ";
				if (coherence) strm2_new +="Cl: ";
				break;
			case 7:
				measures->SetComputeShapeCoefficients( 1 );
				if (!coherence) strm1_new +="Cp: ";
				if (coherence) strm2_new +="Cp: ";
				break;
			case 8:
				measures->SetComputeShapeCoefficients( 2 );
				if (!coherence) strm1_new +="Cs: ";
				if (coherence) strm2_new +="Cs: ";
				break;
			case 9:
				measures->SetComputeDC(0);
				if (!coherence) strm1_new +="D0: ";
				if (coherence) strm2_new +="D0: ";
				break;
			case 10:
				measures->SetComputeDC(1);
				if (!coherence) strm1_new +="D1: ";
				if (coherence) strm2_new +="D1: ";
				break;
			case 11:
				measures->SetComputeDC(2);
				if (!coherence) strm1_new +="D2: ";
				if (coherence) strm2_new +="D2: ";
				break;
			case 12:
				measures->SetComputeDC(3);
				if (!coherence) strm1_new +="D3: ";
				if (coherence) strm2_new +="D3: ";
				break;
			case 13:
				measures->SetComputeDC(4);
				if (!coherence) strm1_new +="D4: ";
				if (coherence) strm2_new +="D4: ";
				break;
			case 14:
				measures->SetComputeDC(5);
				if (!coherence) strm1_new +="D5: ";
				if (coherence) strm2_new +="D5: ";
				break;
			default:
				break;
		}			
		avg_measure.clear();
		FloatImageType::Pointer	measureImage=FloatImageType::New();
			
		typedef itk::CoherenceFilter<FloatImageType, FloatImageType, UCharImageType> coherenceFilterType;
		coherenceFilterType::Pointer coherenceFilter=coherenceFilterType::New();
			
		unsigned long num_voxels;
		for(unsigned i=0; i<tractIndex.size(); ++i){
			measures->SetLabel(tractIndex[i]+1);
				
			measures->Update();
			measureImage=measures->GetOutput();
				
			if(coherence){
				coherenceFilter->SetInput1(measureImage);
				coherenceFilter->SetInput2(measureImage);
				coherenceFilter->SetLabel(tractIndex[i]+1);
				coherenceFilter->SetMask(m_TractRegionsImage);
				coherenceFilter->Update();
				measureImage=coherenceFilter->GetOutput();
			}
			it=itk::ImageRegionConstIterator<FloatImageType>(measureImage, measureImage->GetRequestedRegion() );
				
			avg_measure.push_back(0.0);
			num_voxels=0;
			for( it.GoToBegin(), mit.GoToBegin(); !mit.IsAtEnd(); ++mit, ++it ){
				if(mit.Get()==(tractIndex[i]+1) ){
					avg_measure[i]+=(it.Get());
					num_voxels++;
				}
			}
			
			if(num_voxels>0){
				avg_measure[i]=avg_measure[i]/num_voxels;
			}
			if(!coherence){
				//strm1<<"Region: "<<tractIndex[i]<< " Valor: "<<avg_measure[i]<<std::endl;
				sprintf(sf,"%0.5f ",avg_measure[i]);
				if (strlen(sf) == 8)  strm1_new += " " + std::string(sf);
				if (strlen(sf) == 9)  strm1_new += std::string(sf);
					
			}else{
				//strm2<<"Region: "<<tractIndex[i]<< " Valor: "<<avg_measure[i]<<std::endl;
				sprintf(sf,"%0.5f ",avg_measure[i]);
				if (strlen(sf) == 8)  strm2_new += " " + std::string(sf);
				if (strlen(sf) == 9)  strm2_new += std::string(sf);
			}
			//	std::cout<<"Region "<<tractIndex[i]<<" Medida: "<<avg_measure[i]<<std::endl;
		}
		strm1_new += "\n";
		strm2_new += "\n";
			
	}//End for ind_measure=0...15
	
	if (coherence) {
		m_bufferROIsText->text(strm2_new.c_str());
	} else {
	    m_bufferROIsText->text(strm1_new.c_str());
	}
	
};

void TensorConsole::CorrelateTractRegion(std::vector<unsigned int> tractIndex, unsigned int ind_measure ){
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	
	ComputeMeasuresPointer measuresData=ComputeMeasuresType::New();
	ComputeMeasuresPointer measuresModel=ComputeMeasuresType::New();
	
	measuresData->SetInput((*m_VectorTensorData)[dataId].image);
	measuresData->SetInput2(m_TractRegionsImage);
	
	itk::ImageFileReader<TensorImageType>::Pointer readModel=itk::ImageFileReader<TensorImageType>::New();
	readModel->SetFileName("modelo.vtk");
	itk::ImageFileReader<UCharImageType>::Pointer readRegions=itk::ImageFileReader<UCharImageType>::New();
	readRegions->SetFileName("regiones.vtk");
	
	try{
		readModel->Update();
	}
	catch(itk::ExceptionObject & e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	
	try{
		readRegions->Update();
	}
	catch(itk::ExceptionObject & e){
		fl_alert( e.GetDescription() );
		return;
	}
	
	measuresModel->SetInput(readModel->GetOutput());
	measuresModel->SetInput2(readRegions->GetOutput());
	
	itk::ImageRegionConstIterator<UCharImageType>  mit = itk::ImageRegionConstIterator<UCharImageType>( m_TractRegionsImage, m_TractRegionsImage->GetRequestedRegion() );
	itk::ImageRegionConstIterator<FloatImageType>  it;
	
	ofstream strm("CorrelacionMedidas.txt");
	
	
	std::vector<float> corr_measure;
	for(ind_measure=0; ind_measure<15; ++ind_measure){
		strm<<"Medida "<<ind_measure<<":"<<std::endl;
		switch(ind_measure){
			case 0:
				measuresData->SetComputeFA( );
				measuresModel->SetComputeFA( );
				
				break;
			case 1:
				measuresData->SetComputeRA( );
				measuresModel->SetComputeRA( );				
				break;
			case 2:
				measuresData->SetComputeEigVal( 0 );
				measuresModel->SetComputeRA( );
				break;
			case 3:
				measuresData->SetComputeEigVal( 1 );
				measuresModel->SetComputeEigVal( 1 );
				break;
			case 4:
				measuresData->SetComputeEigVal( 2 );
				measuresModel->SetComputeEigVal( 2 );
				break;
			case 5:
				measuresData->SetComputeMD(  );
				measuresModel->SetComputeMD(  );
				break;
			case 6:
				measuresData->SetComputeDC(0);
				measuresModel->SetComputeDC(0);
				break;
			case 7:
				measuresData->SetComputeDC(1);
				measuresModel->SetComputeDC(1);
				break;
			case 8:
				measuresData->SetComputeDC(2);
				measuresModel->SetComputeDC(2);
				break;
			case 9:
				measuresData->SetComputeDC(3);
				measuresModel->SetComputeDC(3);
				break;
			case 10:
				measuresData->SetComputeDC(4);
				measuresModel->SetComputeDC(4);
				break;
			case 11:
				measuresData->SetComputeDC(5);
				measuresModel->SetComputeDC(5);
				break;
			case 12:
				measuresData->SetComputeShapeCoefficients( 0 );
				measuresModel->SetComputeShapeCoefficients( 0 );
				break;
			case 13:
				measuresData->SetComputeShapeCoefficients( 1 );
				measuresModel->SetComputeShapeCoefficients( 1 );
				break;
			case 14:
				measuresData->SetComputeShapeCoefficients( 2 );
				measuresModel->SetComputeShapeCoefficients( 2 );
				break;
			default:
				measuresData->SetComputeFA();
				measuresModel->SetComputeFA();
				break;
		}			
		corr_measure.clear();
		FloatImageType::Pointer	correlationImage=FloatImageType::New();
		
		/*measureImage->SetRegions(m_TractRegionsImage->GetRequestedRegion());
		 measureImage->SetOrigin(m_TractRegionsImage->GetOrigin());
		 measureImage->SetSpacing(m_TractRegionsImage->GetSpacing());
		 measureImage->Allocate();^*/
		
		TransformType::Pointer transform=TransformType::New();
		transform->SetParameters(m_Parameters);
		
		typedef itk::ResampleImageFilter< FloatImageType, FloatImageType >				 ResampleType;
		typedef itk::NearestNeighborInterpolateImageFunction< FloatImageType, double >	 InterpolatorType;
		
		ResampleType::Pointer			resample		=	ResampleType::New();
		InterpolatorType::Pointer		interpolator    =	InterpolatorType::New();
		resample->SetInterpolator( interpolator );
		resample->SetSize(         m_TractRegionsImage->GetLargestPossibleRegion().GetSize() );
		resample->SetOutputOrigin( m_TractRegionsImage->GetOrigin()                         );
		resample->SetOutputSpacing(m_TractRegionsImage->GetSpacing()		);
		resample->SetTransform(     transform                  );
		
		typedef itk::CorrelationCoefficientFilter<FloatImageType, FloatImageType, UCharImageType> correlationFilterType;
		correlationFilterType::Pointer correlationFilter=correlationFilterType::New();
		/********************************************************************************************/
		//HAY QUE COMPROBAR SI YA SE HA REGISTRADO, Y SI NO, REGISTRAR PARA OBTENER LOS PARAMETROS
		/********************************************************************************************/
		if(!m_HasBeenRegistered){
			transform->SetIdentity();
		}
		
		unsigned long num_voxels;
		for(unsigned i=0; i<tractIndex.size(); ++i){
			measuresData->SetLabel(tractIndex[i]+1);
			measuresModel->SetLabel(tractIndex[i]+1);
			correlationFilter->SetLabel(tractIndex[i]+1);
			try{
				measuresData->Update();
			}catch(itk::ExceptionObject & e){
				fl_alert( e.GetDescription() );
				return;
			}
			try{
				measuresModel->Update();
			}catch(itk::ExceptionObject & e){
				fl_alert( e.GetDescription() );
				return;
			}
			
			resample->SetInput(measuresModel->GetOutput());
			try{
				resample->Update();
			}	catch(itk::ExceptionObject & e){
				fl_alert( e.GetDescription() );
				return;
			}
			
			correlationFilter->SetInput1(measuresData->GetOutput());
			correlationFilter->SetInput2(resample->GetOutput());
			correlationFilter->SetMask(m_TractRegionsImage);
			correlationFilter->Update();
			correlationImage=correlationFilter->GetOutput();
			it=itk::ImageRegionConstIterator<FloatImageType>(correlationImage, correlationImage->GetRequestedRegion() );
			
			corr_measure.push_back(0.0);
			num_voxels=0;
			for( it.GoToBegin(), mit.GoToBegin(); !mit.IsAtEnd(); ++mit, ++it ){
				if(mit.Get()==(tractIndex[i]+1) ){
					corr_measure[i]+=(it.Get());
					num_voxels++;
				}
			}
			if(num_voxels>0){
				corr_measure[i]=corr_measure[i]/num_voxels;
			}
			
			strm<<"Region: "<<tractIndex[i]<< " Valor: "<<corr_measure[i]<<std::endl;
			
			//std::cout<<"Region "<<tractIndex[i]<<" Medida: "<<corr_measure[i]<<std::endl;
		}
	}
	
	
	
};




void TensorConsole::ComputeAngleCoherenceTractRegion(std::vector<unsigned int> tractIndex, unsigned int ind_measure ){
	int dataId = m_tensordataBrowser->value()-1;
	if( dataId<0 || dataId >= (int)m_VectorTensorData->size() ){return;}
	
	typedef itk::AngleCoherenceFilter<TensorImageType, UCharImageType, FloatImageType> coherenceFilterType;
	coherenceFilterType::Pointer coherenceFilter=coherenceFilterType::New();
	
	coherenceFilter->SetInput((*m_VectorTensorData)[dataId].image);
	coherenceFilter->SetMask(m_TractRegionsImage);
	
	itk::ImageRegionConstIterator<UCharImageType>  mit = itk::ImageRegionConstIterator<UCharImageType>( m_TractRegionsImage, m_TractRegionsImage->GetRequestedRegion() );
	itk::ImageRegionConstIterator<FloatImageType >  it;
	
	ofstream strm("Coherencia.txt");
	
	std::vector<float> coh_measure;
	coh_measure.clear();
	
	unsigned long num_voxels;
	for(unsigned i=0; i<tractIndex.size(); ++i){
		coherenceFilter->SetLabel(tractIndex[i]+1);
		try{
			coherenceFilter->Update();
		}catch(itk::ExceptionObject & e){
			fl_alert( e.GetDescription() );
			return;
		}
		
		it=itk::ImageRegionConstIterator<FloatImageType>(coherenceFilter->GetOutput(), coherenceFilter->GetOutput()->GetRequestedRegion() );			
		coh_measure.push_back(0.0);
		num_voxels=0;
		for( it.GoToBegin(), mit.GoToBegin(); !mit.IsAtEnd(); ++mit, ++it ){
			if(mit.Get()==(tractIndex[i]+1) ){
				coh_measure[i]+=(it.Get());
				num_voxels++;
			}
		}
		if(num_voxels>0){
			coh_measure[i]=coh_measure[i]/num_voxels;
		}
		
		strm<<"Region: "<<tractIndex[i]<< " Valor: "<<coh_measure[i]<<std::endl;
		std::cout<<"Region "<<tractIndex[i]<<" Medida: "<<coh_measure[i]<<std::endl;
	}
	
};





void TensorConsole::GetFiberName( unsigned int value, char* nombre ) 
{
	switch (value) {
		case 1:
			sprintf(nombre,"cst left ");break;
		case 2:
			sprintf(nombre,"ml left ");break;
		case 3:
			sprintf(nombre,"scp left ");break;
		case 4:
			sprintf(nombre,"icp left ");break;
		case 5:
			sprintf(nombre,"tap left ");break;
		case 6:
			sprintf(nombre,"atr left ");break;
		case 7:
			sprintf(nombre,"ptr/str left ");break;
		case 8:
			sprintf(nombre,"sfo left ");break;
		case 9:
			sprintf(nombre,"ifo left ");break;
		case 10:
			sprintf(nombre,"slf left ");break;
		case 11:
			sprintf(nombre,"ilf left ");break;
		case 12:
			sprintf(nombre,"fasc unc left ");break;
		case 13:
			sprintf(nombre,"stria terminalis left ");break;
		case 14:
			sprintf(nombre,"cingulum left ");break;
		case 15:
			sprintf(nombre,"cst right ");break;
		case 16:
			sprintf(nombre,"ml right ");break;
		case 17:
			sprintf(nombre,"scp right ");break;
		case 18:
			sprintf(nombre,"icp right ");break;
		case 19:
			sprintf(nombre,"tap right ");break;
		case 20:
			sprintf(nombre,"atr right ");break;
		case 21:
			sprintf(nombre,"ptr/str right ");break;
		case 22:
			sprintf(nombre,"sfo right ");break;
		case 23:
			sprintf(nombre,"ifo right ");break;
		case 24:
			sprintf(nombre,"slf right ");break;
		case 25:
			sprintf(nombre,"ilf right ");break;
		case 26:
			sprintf(nombre,"fasc unc right ");break;
		case 27:
			sprintf(nombre,"stria terminalis ");break;
		case 28:
			sprintf(nombre,"cingulum right ");break;
		case 29:
			sprintf(nombre,"corpus callosum ");break;
		case 30:
			sprintf(nombre,"fornix ");break;
		case 31:
			sprintf(nombre,"mcp ");break;
		default:
			break;
			
	}
}

void TensorConsole::GetColorFromSeed( unsigned int seedvalue ) 
{
	
	switch (seedvalue ) {
		case 1:
			m_color_r->value(1);m_color_g->value(0);m_color_b->value(0);break;
		case 2:
			m_color_r->value(1);m_color_g->value(1);m_color_b->value(0);break;
		case 3:
			m_color_r->value(0.9);m_color_g->value(0);m_color_b->value(0.5);break;
		case 4:
			m_color_r->value(0.2);m_color_g->value(0.2);m_color_b->value(0.8);break;
		case 5:
			m_color_r->value(0.1);m_color_g->value(0.5);m_color_b->value(0.1);break;
		case 6:
			m_color_r->value(0);m_color_g->value(1);m_color_b->value(0);break;
		case 7:
			m_color_r->value(1);m_color_g->value(1);m_color_b->value(0.6);break;
		case 8:
			m_color_r->value(0.9);m_color_g->value(0.3);m_color_b->value(0.9);break;
		case 9:
			m_color_r->value(0.7);m_color_g->value(0.3);m_color_b->value(0.2);break;
		case 10:
			m_color_r->value(0.9);m_color_g->value(0.9);m_color_b->value(0.9);break;
		case 11:
			m_color_r->value(0.5);m_color_g->value(1);m_color_b->value(0.5);break;
		case 12:
			m_color_r->value(0.1);m_color_g->value(0.4);m_color_b->value(0.4);break;
		case 13:
			m_color_r->value(0.6);m_color_g->value(0.6);m_color_b->value(0.1);break;
		case 14:
			m_color_r->value(0);m_color_g->value(0);m_color_b->value(1);break;
		case 15:
			m_color_r->value(0.6);m_color_g->value(0.5);m_color_b->value(1);break;
		case 16:
			m_color_r->value(0.6);m_color_g->value(0.1);m_color_b->value(0.1);break;
		case 17:
			m_color_r->value(0.5);m_color_g->value(0);m_color_b->value(0);break;
		case 18:
			m_color_r->value(0.3);m_color_g->value(0.7);m_color_b->value(0.3);break;
		case 19:
			m_color_r->value(0);m_color_g->value(0);m_color_b->value(1);break;
		case 20:
			m_color_r->value(1);m_color_g->value(1);m_color_b->value(1);break;
		case 21:
			m_color_r->value(0.8);m_color_g->value(0.2);m_color_b->value(0.7);break;
		case 22:
			m_color_r->value(0.7);m_color_g->value(0.5);m_color_b->value(0.2);break;
		case 23:
			m_color_r->value(0.2);m_color_g->value(0.7);m_color_b->value(0.5);break;
		case 24:
			m_color_r->value(0.5);m_color_g->value(1);m_color_b->value(1);break;
		case 25:
			m_color_r->value(1);m_color_g->value(0.5);m_color_b->value(1);break;
		case 26:
			m_color_r->value(0.5);m_color_g->value(0);m_color_b->value(1);break;
		case 27:
			m_color_r->value(1);m_color_g->value(0.5);m_color_b->value(0);break;
		case 28:
			m_color_r->value(1);m_color_g->value(0);m_color_b->value(0.5);break;
		case 29:
			m_color_r->value(0);m_color_g->value(0.5);m_color_b->value(1);break;
		case 30:
			m_color_r->value(0.2);m_color_g->value(1);m_color_b->value(0.6);break;
		case 31:
			m_color_r->value(0.6);m_color_g->value(0.2);m_color_b->value(0.2);break;
		default:
			break;
	}
	
}

void TensorConsole::RegisterSeedsToData(void){
	
	int dataId = m_tensordataBrowser->value()-1;
	if(dataId<0){return;}
	
	typedef itk::ImageFileReader< TensorImageType  >			ImageReaderType;
	ImageReaderType::Pointer		reader		=		ImageReaderType::New();
	
	typedef itk::ImageFileReader< UCharImageType >			LabelImageReaderType; 
	LabelImageReaderType::Pointer labelReader =				LabelImageReaderType::New();
	
	reader->SetFileName("modelo.vtk" );				// ¿Se deja elegir el modelo o usamos siempre este?
	labelReader->SetFileName( "ROIs.mhd");		// Esto tiene que cambiarse cuando se tenga la imagen segmentada
	
	TensorImageType::Pointer tensor = (*m_VectorTensorData)[dataId].image;
	
	m_register	=	FaTensorRegistrationType::New();
	m_register->SetFixed(tensor);
	
	try{
		reader->Update( );
	}
	catch ( itk::ExceptionObject & e )
	{
		fl_alert( e.GetDescription() );
		return;
	}
	m_register->SetMoving(reader->GetOutput());
	
	try{
		labelReader->Update();
	}
	catch ( itk::ExceptionObject & e )
	{
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_register->SetLabels(labelReader->GetOutput());
	
	try{
		m_register->Start();
	}
	catch ( itk::ExceptionObject & e )
	{
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_SeedRegionsImage=m_register->GetOutput();
	m_Parameters=m_register->GetParameters();
	
	ImageViewer[m_activeinput]->SetInputOverlay(m_SeedRegionsImage);
	return;
	
};

void TensorConsole::LoadROIS(void){
	
	typedef itk::ImageFileReader< UCharImageType >			LabelImageReaderType; 
	LabelImageReaderType::Pointer labelReader =				LabelImageReaderType::New();
	labelReader->SetFileName( "ROIs.mhd");		// Esto tiene que cambiarse cuando se tenga la imagen segmentada
	
	try{
		labelReader->Update();
	}
	catch ( itk::ExceptionObject & e )
	{
		fl_alert( e.GetDescription() );
		return;
	}
	
	m_SeedRegionsImage=labelReader->GetOutput();
	ImageViewer[m_activeinput]->SetInputOverlay(m_SeedRegionsImage);
	m_imageOverlay=m_SeedRegionsImage;
	m_HasBeenRegistered = true;
}

void TensorConsole::FilterDWI( DWImagesType::Pointer &dwimage) {  
	
	m_dwifilter->SetInput( dwimage );
	
	//==========================================================================
	// Set the parameters to the filter:
	DWIFilterType::InputSizeType radius;
	//unsigned int rade = 1;
	radius.Fill( (unsigned int)m_rade->value() );
	//-----------------------
	//radius[2] = 0;
	//-----------------------
	m_dwifilter->SetRadiusEstimation( radius );
	//unsigned int radf = 2;
	radius.Fill( (unsigned int)m_radf->value() );
	radius[2] = 0;
	m_dwifilter->SetRadiusFiltering( radius );
	//unsigned int iter1 = 1;
	
	m_dwifilter->SetIterations( (unsigned int)m_iter1->value() );
	m_dwifilter->SetUseAbsoluteValue( false );
	m_dwifilter->SetKeepValue( true );
	m_dwifilter->SetMinimumNumberOfUsedVoxelsFiltering( 5 );
	m_dwifilter->SetMinimumNoiseSTD( 0 );
	m_dwifilter->SetMaximumNoiseSTD( 32000 );
	m_dwifilter->SetFirstBaseline( 0 );
	m_dwifilter->SetChannels( dwimage->GetVectorLength() );
	//==========================================================================
	
	try
	{
		if( (unsigned int)m_iter1->value()>0 ){
			std::cerr << "Filtering..." << std::endl;
			m_dwifilter->Update();
			std::cerr << "Done!" << std::endl;
			m_dwifilter->GetOutput()->SetNumImages( dwimage->GetNumImages() );
			m_dwifilter->GetOutput()->SetNumDWImages( dwimage->GetNumDWImages() );
			m_dwifilter->GetOutput()->SetNumBaselines( dwimage->GetNumBaselines() );
			m_dwifilter->GetOutput()->SetIndexesDWImages( dwimage->GetIndexesDWImages() );
			m_dwifilter->GetOutput()->SetIndexesBaselines( dwimage->GetIndexesBaselines() );
			m_dwifilter->GetOutput()->SetBValues( dwimage->GetBValues() );
			DWImagesType::DiffusionDirectionsType directions;
			dwimage->GetDiffusionDirections( directions );
			m_dwifilter->GetOutput()->SetDiffusionDirections( directions );
		}
	}
	catch ( itk::ExceptionObject & e )
	{
		fl_alert( e.GetDescription() );
		return;
	}
	dwimage = m_dwifilter->GetOutput();
}

void TensorConsole::EstimateTensor( DWImagesType::Pointer dwimage, TensorImageType::Pointer &tensor, InputImageType::Pointer &t2image, bool mask) {  
	
	m_dtiEstimator->SetInput( dwimage );
	// The number of iterations; by default, only one iteration (simple WLS) is
	// performed:
	unsigned int iter2 = 1;
	
	m_dtiEstimator->SetIterations( iter2 );
	std::cerr << "Configuring DTI estimator..." << std::endl;
	m_dtiEstimator->Configure();
	std::cerr << "Done." << std::endl;
	std::cerr << "Computing the tensor..." << std::endl;
	
	m_dtiEstimator->SetComputeT2(true);
	
	if (mask) {
	  try{
		m_dtiEstimator->ComputeMask();
	  }
	  catch ( itk::ExceptionObject & e )
	  {
		fl_alert( e.GetDescription() );
		return;
	  }
	}
	
	try
	{
		m_dtiEstimator->Update();
		std::cout << "Done" << std::endl;
	}
	catch ( itk::ExceptionObject & e )
	{
		fl_alert( e.GetDescription() );
		return;
	}
	tensor  = m_dtiEstimator->GetOutput();
	t2image = m_dtiEstimator->GetT2();
}

void TensorConsole::GetComponentFromDWI( DWImagesType::Pointer dwimage, unsigned int n) {  
	
	int id = m_DWIdataBrowser->value()-1;
	if (id<0){return;}
	
	InputImageType::Pointer image	= InputImageType::New(); 
	image->SetRegions((*m_VectorDWIData)[id].image->GetLargestPossibleRegion());
	image->SetOrigin((*m_VectorDWIData)[id].image->GetOrigin());
	image->SetSpacing((*m_VectorDWIData)[id].image->GetSpacing());
	image->Allocate();
	
	(*m_VectorDWIData)[id].image->GetImageComponent(n, image);
	
	try
	{
		(*m_VectorDWIData)[id].image->Update();
		std::cout << "Done" << std::endl;
	}
	catch ( itk::ExceptionObject & e )
	{
		fl_alert( e.GetDescription() );
		return;
	}
	
	char nombre[200];
	sprintf( nombre, "Component_%02d",n);
	
	m_VectorData->generateNewData( image, nombre );
	ImageViewer[m_activeinput]->SetImage( image );
	image = NULL;
	return;
}

void TensorConsole::Mascara( TensorImageType::Pointer tensor, InputImageType::Pointer mask) {  
	
	ConstIteratorType imageIt(mask,mask->GetLargestPossibleRegion());
	TensorIteratorType tensorIt(tensor,tensor->GetLargestPossibleRegion());
	for (tensorIt.GoToBegin(),imageIt.GoToBegin();!tensorIt.IsAtEnd();++tensorIt, ++imageIt) {
		tensorIt.Set(tensorIt.Get()*imageIt.Get());
	}
	
}

void TensorConsole::ApplyMask(  TensorImageType::Pointer &tensor ) {  
	int dataId = m_dataBrowser->value()-1;
	this->Mascara(tensor, (*m_VectorData)[dataId].image);
	std::cout << "Calculada mascara!" << std::endl;
}

void TensorConsole::verGlifosTract() {

	if (!mostrarGlifos->value()) return;

	int dataId = m_tensordataBrowser->value()-1;

	vtkTensorGlyphDTI *gen = new vtkTensorGlyphDTI();
	gen->SetInput((*m_VectorTensorData)[dataId].image);
	gen->SetInputPoints(m_puntosTract);
	gen->SetBounds(0,0,0,0,0,0);
	gen->SetScaleFactor(valorEscala->value());
	gen->SetPhiResolution(valorPhiResolution->value());
	gen->SetThetaResolution(valorThetaResolution->value());
	gen->SetGamma(valorGamma->value());

	if (verCuboides->value()) gen->SetGlyphType(vtkTensorGlyphDTI::CUBOID);
	else if (verSupercuadricas->value()) gen->SetGlyphType(vtkTensorGlyphDTI::SUPERQUADRIC);
	else gen->SetGlyphType(vtkTensorGlyphDTI::ELLIPSOID);

	if (colorRA->value()) gen->SetColorMode(vtkTensorGlyphDTI::COLOR_BY_RA);
	else if (colorFA->value()) gen->SetColorMode(vtkTensorGlyphDTI::COLOR_BY_FA);
	else if (colorCl->value()) gen->SetColorMode(vtkTensorGlyphDTI::COLOR_BY_CL);

	borrarGlifos(3);

	vtkPolyData *glifos = vtkPolyData::New();

	gen->GetOutput(glifos);

	m_tractActor = vtkActor::New();

	ImageViewer3D->ConnectMapper(glifos, m_tractActor);
	ImageViewer3D->SetScalarRange(m_tractActor, 0, 1);	

	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();

	glifos->Delete();
}

void TensorConsole::borrarGlifosStrain() {

	if (m_activeActorStrain) {
		ImageViewer3D->HideModel(m_activeActorStrain);
		m_activeActorStrain->GetMapper()->Delete();
		m_activeActorStrain->Delete();
		m_activeActorStrain=NULL;
	}

}

void TensorConsole::verGlifosStrain() {

	verGlifosStrain(ImageViewerStrain3D->planeWidgetZ->GetSliceIndex(), this->tiempo);

}

void TensorConsole::verGlifosStrain(int planoZ, int tiempo) {

	this->tiempo = tiempo;

	borrarGlifosStrain();

	int dataId = m_strainTensorDataBrowser->value()-1;
	if( dataId<0){return;}

	STImageType::Pointer imagen = (*m_VectorSTData)[dataId].image;
	DeformImageType::Pointer deform = (*m_VectorSTData)[dataId].deform_image;

	vtkPoints *puntos = vtkPoints::New();
	puntos->Allocate(0);	

	STImageType::IndexType indice;
	indice[2] = planoZ;
	indice[3] = tiempo;

	int xSize = imagen->GetRequestedRegion().GetSize()[0];
	int ySize = imagen->GetRequestedRegion().GetSize()[1];

	STImageType::PointType origin = imagen->GetOrigin();
	STImageType::SpacingType spacing = imagen->GetSpacing();

	float escala = 1.0;

	// Si se aumenta la resolución de glifos, se calculan los puntos adicionales donde mostrarlos
	if(strainMuestreo2->value()) {

		escala = 0.5;

		for (int i=0; i<xSize; i++) {

			indice[0] = i;

			for (int j=0; j<ySize; j++) {
		
				indice[1] = j;

				if (i<xSize-1) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.5),origin[1]+spacing[1]*(j),origin[2]+spacing[2]*planoZ);
				}
				if (j<ySize-1) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i),origin[1]+spacing[1]*(j+0.5),spacing[2]*planoZ);
				}
				if ( (i<xSize-1) && (j<ySize-1) ) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.5),origin[1]+spacing[1]*(j+0.5),spacing[2]*planoZ);
				}

			}
		}
	}

	if(strainMuestreo3->value()) {
	
		escala = 0.33;

		for (int i=0; i<xSize; i++) {

			indice[0] = i;

			for (int j=0; j<ySize; j++) {
		
				indice[1] = j;

				if (i<xSize-1) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.33),origin[1]+spacing[1]*(j),origin[2]+spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.66),origin[1]+spacing[1]*(j),origin[2]+spacing[2]*planoZ);
				}
				if (j<ySize-1) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i),origin[1]+spacing[1]*(j+0.33),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i),origin[1]+spacing[1]*(j+0.66),spacing[2]*planoZ);
				}
				if ( (i<xSize-1) && (j<ySize-1) ) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.33),origin[1]+spacing[1]*(j+0.33),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.33),origin[1]+spacing[1]*(j+0.66),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.66),origin[1]+spacing[1]*(j+0.33),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.66),origin[1]+spacing[1]*(j+0.66),spacing[2]*planoZ);
				}

			}
		}
	}

	if(strainMuestreo4->value()) {

		escala = 0.25;

		for (int i=0; i<xSize; i++) {

			indice[0] = i;

			for (int j=0; j<ySize; j++) {
		
				indice[1] = j;

				if (i<xSize-1) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.25),origin[1]+spacing[1]*(j),origin[2]+spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.5),origin[1]+spacing[1]*(j),origin[2]+spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.75),origin[1]+spacing[1]*(j),origin[2]+spacing[2]*planoZ);
				}
				if (j<ySize-1) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i),origin[1]+spacing[1]*(j+0.25),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i),origin[1]+spacing[1]*(j+0.5),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i),origin[1]+spacing[1]*(j+0.75),spacing[2]*planoZ);
				}
				if ( (i<xSize-1) && (j<ySize-1) ) {
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.25),origin[1]+spacing[1]*(j+0.25),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.25),origin[1]+spacing[1]*(j+0.5),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.25),origin[1]+spacing[1]*(j+0.75),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.5),origin[1]+spacing[1]*(j+0.25),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.5),origin[1]+spacing[1]*(j+0.5),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.5),origin[1]+spacing[1]*(j+0.75),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.75),origin[1]+spacing[1]*(j+0.25),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.75),origin[1]+spacing[1]*(j+0.5),spacing[2]*planoZ);
					puntos->InsertNextPoint(origin[0]+spacing[0]*(i+0.75),origin[1]+spacing[1]*(j+0.75),spacing[2]*planoZ);
				}

			}
		}
	}

	vtkTensorGlyphStrain *glifos = new vtkTensorGlyphStrain();
	glifos->SetInput(imagen);
	glifos->SetDeformImage(deform);
	glifos->SetPlanoZ(planoZ);
	glifos->SetTiempo(tiempo);
	glifos->SetInputPoints(puntos);

	glifos->SetColorMode(vtkTensorGlyphStrain::DEF);

	if (strainST0->value())
		glifos->SetColorMode(vtkTensorGlyphStrain::ST0);
	else if (strainST1->value())
		glifos->SetColorMode(vtkTensorGlyphStrain::ST1);
	else if (strainST2->value())
		glifos->SetColorMode(vtkTensorGlyphStrain::ST2);
	else if (strainEIG0->value())
		glifos->SetColorMode(vtkTensorGlyphStrain::EIG0);
	else if (strainEIG1->value())
		glifos->SetColorMode(vtkTensorGlyphStrain::EIG1);
	else glifos->SetColorMode(vtkTensorGlyphStrain::DEF);

	glifos->SetScaleFactor(escala * strainEscala->value());

	m_activeActorStrain = vtkActor::New();
	ImageViewerStrain3D->ConnectMapper(glifos->GetOutput(), m_activeActorStrain);
	ImageViewerStrain3D->SetScalarRange(m_activeActorStrain, rangeStrainMin, rangeStrainMax);

	if (!m_scalarBarStrain) {
		m_scalarBarStrain = vtkScalarBarActor::New();
		m_scalarBarStrain->SetHeight(0.3);
		m_scalarBarStrain->SetWidth(0.06);
	}

	vtkLookupTable *lut = vtkLookupTable::New();
	lut->SetRange(rangeStrainMin,rangeStrainMax);
	lut->Build();

	m_scalarBarStrain->SetLookupTable(lut);
	m_activeActorStrain->GetMapper()->SetLookupTable(lut);

	ImageViewerStrain3D->AddViewProp(m_scalarBarStrain);


}

void TensorConsole::SetScalarRangeStrain(double min, double max) {
	rangeStrainMin = min;
	rangeStrainMax = max;
	ImageViewerStrain3D->SetScalarRange(m_activeActorStrain, min, max);
}

void TensorConsole::actualizarGlifosDTI() {

	if (m_planoActivoX) verGlifos(0);
	if (m_planoActivoY) verGlifos(1);
	if (m_planoActivoZ) verGlifos(2);
	if (m_tractActiva) verGlifosTract();

	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();
}

void TensorConsole::imagenActiva(int numPlano,bool activo) {

	switch (numPlano) {
		case 0: m_planoActivoX = activo;
			break;

		case 1: m_planoActivoY = activo;
			break;

		case 2: m_planoActivoZ = activo;
			break;

		case 3: m_tractActiva = activo;
			break;

	}

}

void TensorConsole::glifosActivos() {

	if (m_planoActivoX) verGlifos(0);
	
	if (m_planoActivoY) verGlifos(1);

	if (m_planoActivoZ) verGlifos(2);

	if (m_tractActiva) verGlifosTract();

	vtkLookupTable *lut = vtkLookupTable::New();
	lut->Build();

	m_scalarBar = vtkScalarBarActor::New();
	m_scalarBar->SetLookupTable(lut);
	m_scalarBar->SetHeight(0.25);
	m_scalarBar->SetWidth(0.05);

	if (colorRA->value())
		m_scalarBar->SetTitle("RA");

	else if (colorFA->value())				
		m_scalarBar->SetTitle("FA");

	else if (colorCl->value())
		m_scalarBar->SetTitle("Cl");

	ImageViewer3D->AddViewProp(m_scalarBar);


}

void TensorConsole::cambiarOpacidad(float value) {
	
	if (m_activeActorX)
		ImageViewer3D->ChangeOpacity(m_activeActorX, value);

	if (m_activeActorY)
		ImageViewer3D->ChangeOpacity(m_activeActorY, value);

	if (m_activeActorZ)
		ImageViewer3D->ChangeOpacity(m_activeActorZ, value);

	if (m_tractActor)
		ImageViewer3D->ChangeOpacity(m_tractActor, value);


}

void TensorConsole::borrarGlifos(int numImagen) {
	
	switch (numImagen) {
	
		case 0:	if (m_activeActorX) {
				ImageViewer3D->HideModel(m_activeActorX);
				m_activeActorX->GetMapper()->Delete();
				m_activeActorX->Delete();
				m_activeActorX=NULL;
			}
			break;

		case 1:	if (m_activeActorY) {
				ImageViewer3D->HideModel(m_activeActorY);
				m_activeActorY->GetMapper()->Delete();
				m_activeActorY->Delete();
				m_activeActorY=NULL;
			}
			break;

	
		case 2:	if (m_activeActorZ) {
				ImageViewer3D->HideModel(m_activeActorZ);
				m_activeActorZ->GetMapper()->Delete();
				m_activeActorZ->Delete();
				m_activeActorZ=NULL;
			}
			break;

		case 3:	if (m_tractActor) {
				ImageViewer3D->HideModel(m_tractActor);
				m_tractActor->GetMapper()->Delete();
				m_tractActor->Delete();
				m_tractActor=NULL;
			}
			break;

		case 4:	if (m_scalarBar) {
				ImageViewer3D->RemovePropA(m_scalarBar);
				m_scalarBar->Delete();
				m_scalarBar=NULL;
			}
			break;

	}

}

void TensorConsole::verGlifos(int numImagen) {

	if (!mostrarGlifos->value()) return;

	int dataId = m_tensordataBrowser->value()-1;

	int i, j, k, offset, jOffset, kOffset, iMax, jMax, kMax;
	float sp, s, x, y, z;
	double *m[3], w[3], *v[3];
	double m0[3], m1[3], m2[3];
	double v0[3], v1[3], v2[3];
	double factor;

	int xMin, xMax, yMin, yMax, zMin, zMax;

	xMin = xMinimo->value();
	xMax = xMaximo->value();
	yMin = yMinimo->value();
	yMax = yMaximo->value();
	zMin = zMinimo->value();
	zMax = zMaximo->value();

	vtkTensorGlyphDTI *gen = new vtkTensorGlyphDTI();
	gen->SetInput((*m_VectorTensorData)[dataId].image);
	gen->SetScaleFactor(valorEscala->value());
	gen->SetPhiResolution(valorPhiResolution->value());
	gen->SetThetaResolution(valorThetaResolution->value());
	gen->SetGamma(valorGamma->value());

	if (verCuboides->value()) gen->SetGlyphType(vtkTensorGlyphDTI::CUBOID);
	else if (verSupercuadricas->value()) gen->SetGlyphType(vtkTensorGlyphDTI::SUPERQUADRIC);
	else gen->SetGlyphType(vtkTensorGlyphDTI::ELLIPSOID);

	if (colorRA->value()) gen->SetColorMode(vtkTensorGlyphDTI::COLOR_BY_RA);
	else if (colorFA->value()) gen->SetColorMode(vtkTensorGlyphDTI::COLOR_BY_FA);
	else if (colorCl->value()) gen->SetColorMode(vtkTensorGlyphDTI::COLOR_BY_CL);

	if (filtrarFA->value()) gen->SetFilterMode(vtkTensorGlyphDTI::FILTER_BY_FA);
	else if (filtrarCl->value()) gen->SetFilterMode(vtkTensorGlyphDTI::FILTER_BY_CL);
	else if (filtrarCs->value()) gen->SetFilterMode(vtkTensorGlyphDTI::FILTER_BY_CS);

	gen->SetFilterThreshold(valorFiltro->value());

	vtkPolyData *glifos = vtkPolyData::New();

	switch (numImagen) {
		
		case 0: borrarGlifos(0);
			i=ImageViewer3D->planeWidgetX->GetSliceIndex();
			gen->SetBounds(i,i, yMin,yMax, zMin,zMax);
			gen->GetOutput(glifos);
			m_activeActorX = vtkActor::New();
			ImageViewer3D->ConnectMapper(glifos, m_activeActorX);
			ImageViewer3D->SetScalarRange(m_activeActorX, 0, 1);
			break;

		case 1: borrarGlifos(1);
			j=ImageViewer3D->planeWidgetY->GetSliceIndex();
			gen->SetBounds(xMin,xMax, j,j, zMin,zMax);
			gen->GetOutput(glifos);
			m_activeActorY = vtkActor::New();
			ImageViewer3D->ConnectMapper(glifos, m_activeActorY);
			ImageViewer3D->SetScalarRange(m_activeActorY, 0, 1);
			break;

		case 2: borrarGlifos(2);
			k=ImageViewer3D->planeWidgetZ->GetSliceIndex();
			gen->SetBounds(xMin,xMax, yMin,yMax, k,k);
			gen->GetOutput(glifos);
			m_activeActorZ = vtkActor::New();
			ImageViewer3D->ConnectMapper(glifos, m_activeActorZ);
			ImageViewer3D->SetScalarRange(m_activeActorZ, 0, 1);
			break;

	}

	Fl::check();
	ImageViewer3D->redraw();
	Fl::check();

	glifos->Delete();

}


