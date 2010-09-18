/*=========================================================================

  Program:   TensorConsole.h
  Language:  C++
  Date:      7-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#ifndef _TensorConsole_h
#define _TensorConsole_h

#include "TensorGUI.h"
#include "Viewer3D.h"
#include "geodesicPath3D.h"
#include "tensor/itkTractography.h"
#include "tensor/itkLMMSEVectorImageFilter.h"
#include "tensor/itkDTIEstimateTensorFilter.h"
#include "tensor/itkDTITensor.h"
#include "tensor/itkDWImages.h"
#include "tensor/itkComputeDTIScalars.h"
#include "tensor/itkComputeMeasuresInRegions.h"
#include "tensor/itkFaDTIRegistration.h"
#include "tensor/itkCorrelationCoefficientFilter.h"
#include "tensor/itkCoherenceFilter.h"
#include "tensor/itkAngleCoherenceFilter.h"

#include <itkFastMarchingImageFilter.h>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Text_Display.H>
#include "MyfltkImageViewer/MyfltkImageViewer.h"
#include "MyfltkImageViewer/MyfltkColorImageViewer.h"
#include <itkRGBPixel.h>
// Usando clase de tensores de ITK:
//#include <itkDiffusionTensor3D.h>
#include <vtkPolyData.h>
#include "VolumesContainer.h"
#include "itkImageFileReader.h"

#include <vtkActor.h>
#include <vtkScalarBarActor.h>
#include <vtkDataSet.h>
#include <vtkLODActor.h>

#include <strain/itkStrainTensor.h>
#include <strain/itkComputeStrainScalars.h>
#include <strain/itkComputeDeformFilter.h>
#include <itkVTKImageExport.h>
#include <itkExtractImageFilter.h>

class TensorConsole : public TensorGUI
{
public:
  enum { Dimension =  3 };
  typedef float                                      InputPixelType;
  typedef float                                      FloatPixelType;
  typedef itk::Image< FloatPixelType, Dimension >    FloatImageType;
  typedef itk::Image< FloatPixelType, 2 >            FloatImageType2D;
  //typedef itk::DWImages< InputPixelType, 3 >         DWImagesType;
  typedef itk::FastMarchingImageFilter<FloatImageType, FloatImageType >   FastMarchingFilter3DType;
  typedef itk::DTITensor<float>								TensorPixelType;

  typedef itk::FixedArray<float,2>	DeformPixelType;
  typedef itk::Image<DeformPixelType,4>	DeformImageType;
  typedef itk::Image<itk::FixedArray<float,2>, 3>	DeformImage3DType;

  // Usando clase de tensores de ITK:
  //typedef itk::DiffusionTensor3D<float> TensorPixelType;
  typedef InputImageType::PointType                           InputPointType;
  typedef itk::Image<TensorPixelType, Dimension>              TensorImageType;
  typedef itk::ImageRegionIterator<InputImageType>            IteratorType;
  typedef itk::RGBPixel<float>                                RGBPixelType;
  typedef itk::Image< RGBPixelType, 3 >                       ImageRGBType;
  typedef itk::Image< unsigned char, 3 >   	                  UCharImageType;
  typedef itk::ImageRegionIterator<ImageRGBType>              RGBIteratorType;
  typedef itk::ImageRegionConstIterator<UCharImageType>       UCharIteratorType;
  typedef itk::ImageRegionConstIterator<TensorImageType>      ConstTensorIteratorType;
  typedef itk::ImageRegionConstIterator<InputImageType>       ConstIteratorType;
  typedef itk::ImageRegionIterator<TensorImageType>           TensorIteratorType;
  typedef itk::ImageFileWriter< FloatImageType2D >            FloatWriter2DType;

  typedef FastMarchingFilter3DType::NodeType                  NodeType;
  typedef FastMarchingFilter3DType::NodeContainer             NodeContainer;
  typedef fltk::MyColorImageViewer<float,float>               MyColorViewerType;
	
  typedef itk::geodesicPath3D< InputImageType, InputImageType >      geodesicPath3DType;

  typedef itk::VectorImage<unsigned long> VectorImageType;
  typedef itk::Tractography< TensorImageType, VectorImageType >      TractographyType;
  typedef itk::LMMSEVectorImageFilter<DWImagesType,DWImagesType>     DWIFilterType;
  typedef itk::DTIEstimateTensorFilter<DWImagesType,TensorImageType> DTIEstimatorType;
	
  typedef itk::ComputeDTIScalars<TensorImageType,InputImageType>     ComputeScalarsType;
  typedef ComputeScalarsType::Pointer                                ComputeScalarsPointer;

  typedef itk::ComputeMeasuresInRegions<TensorImageType, UCharImageType, FloatImageType>     ComputeMeasuresType;
  typedef ComputeMeasuresType::Pointer					ComputeMeasuresPointer;
  
  typedef itk::FaDTIRegistration<TensorImageType,TensorImageType,UCharImageType>		FaTensorRegistrationType;
  
  typedef itk::D3BSplineMaskTransform< double, 3>									TransformType;
  typedef TransformType::ParametersType												ParametersType;

  enum { StrainDimension =  4 };
  typedef itk::StrainTensor<float> 		STPixelType;
  typedef itk::Image<STPixelType, StrainDimension> 	STImageType;
  typedef itk::Image< STPixelType, 3 >    				STImageType3D;
  typedef itk::ComputeStrainScalars <STImageType3D,InputImageType>     	ComputeStrainScalarsType;
  typedef ComputeStrainScalarsType::Pointer                            	ComputeStrainScalarsPointer;
  typedef itk::ComputeDeformFilter <DeformImage3DType,InputImageType>	ComputeDeformFilterType;
  typedef ComputeDeformFilterType::Pointer				ComputeDeformFilterPointer;
  typedef itk::VTKImageExport< InputImageType> VTKexportType;
  VTKexportType::Pointer       m_VTKexporter; 
  typedef itk::ExtractImageFilter< DeformImageType, DeformImage3DType >		ExtractFilterType;
  typedef enum StrainScalars{INV,EIG0,EIG1,ST0,ST1,ST2} StrainScalars;


  // Funciones para tensores
  TensorConsole(int X, int Y, int W, int H, const char *L = 0);
  virtual ~TensorConsole();

  //void SetActiveInput(unsigned int &activeinput);
  void SetModelDataBrowser(Fl_Browser *modeldataBrowserIn);
  void SetDWIDataBrowser(Fl_Browser *DWIdataBrowserIn);
  void SetStrainTensorDataBrowser(Fl_Browser *STDataBrowserIn);
  //void SetImageViewer3D(Viewer3D* ImageViewer3DIn);
  void SetImageColorViewer(MyColorViewerType* ImageColorViewer3DIn);
  void SetImageOverlay(UCharImageType* m_imageOverlayIn);
  void SetSliders(Fl_Value_Slider **sliders); 
  void SetPanelMessage(Fl_Double_Window *panelmessage); 
  void SetOrientationText(Fl_Text_Display *orientationTextIn);
  void SetOriginText(Fl_Text_Display *originTextIn);
  void SetMinScalarRange(Fl_Value_Slider *min_scalar_range_In);
  void SetMaxScalarRange(Fl_Value_Slider *max_scalar_range_In);
  void SetActiveStrainImage(vtkImageData *image);
		
	void InfoImagen();
	void InitializePoints();
	void InitializePoints(unsigned int);
	
  //void RungeKuttaTractographyWithConnection();
  void RungeKuttaTractography(void); 
  void BruteForceTractography();
  void SelectFiberTracts(int); 
  void SelectAutoFiberTracts(int label1);
  void SelectFiberTractsConnection(int, int, std::vector<int>);

  void geodesicPath3D(void);
  void CreateStreamLine(void);
  void CreateStreamLineFloat(void);
  void CreateStreamLineStraight(void);
  void CreateStreamLineNew(void);
  void ComputeScalarValue(  unsigned int scalar, const char* nameScalar, unsigned int window );
  void ColorOrientation( void );
  
  
  void SetScalarsSize(vtkPolyData* streamLineList);
  void SetScalarsDistance(vtkPolyData* streamLineList);
  void SetScalarsMajorEig(vtkPolyData* streamLineList);
  void SetScalarsFA(vtkPolyData* streamLineList);
  void SetScalarsMD(vtkPolyData* streamLineList);
  void SetScalarsCurvature(vtkPolyData* streamLineList);
  void RetainNearbyFibers(vtkPolyData* streamLineList, float factor, int labelROI);
  void OrderFibersCoronal(vtkPolyData* streamLineList, int z);
  void OrderFibersCoronal2(vtkPolyData* streamLineList, int y);
  void CutFibers(vtkPolyData* streamLineList, int orientation, int value, int direccion);
  void ResizeLargeFibers(vtkPolyData* streamLineList, float factor);
  void RemoveLargeFibers(vtkPolyData* streamLineList, float factor);	
  void ComputeFiberMeasures(bool save_stats_file);
  void ComputeFibersSize(vtkPolyData* streamLineList, bool save_stats_file);	
  void ComputeFibersAngle(vtkPolyData* streamLineList);
  void ComputeFibersFA(vtkPolyData* streamLineList, bool save_profile, bool save_stats_file);
  void ComputeFibersRA(vtkPolyData* streamLineList, bool save_stats_file);
  void ComputeFibersMD(vtkPolyData* streamLineList, bool save_stats_file);
  void ComputeFibersEigen(vtkPolyData* streamLineList, bool save_stats_file);	
  void ComputeFibersShapeCoefficients(vtkPolyData* streamLineList, bool save_stats_file);
  void ComputeFibersTensorComponents(vtkPolyData* streamLineList, bool save_stats_file);
  void ComputeFibersFAandProject(vtkPolyData* streamLineList);
  void SetFibersColor(int value);

  void CreateAutomaticFiberTract(std::vector<unsigned int>);
  void CreateAutomaticFiberTractRuben(std::vector<unsigned int>);
  void ProjectROIsFromFibers(vtkPolyData* streamLineList);
  void EvaluateTractRegion(std::vector<unsigned int>, std::vector<unsigned int>, bool coherence);
  void CorrelateTractRegion(std::vector<unsigned int>, unsigned int);
  void ComputeAngleCoherenceTractRegion(std::vector<unsigned int>, unsigned int);
  
  void GetFiberName( unsigned int value, char* nombre);
  void GetColorFromSeed( unsigned int seedvalue );
  void Mascara( TensorImageType::Pointer tensor, InputImageType::Pointer mask);
  void RegisterSeedsToData(void); 
  void FilterDWI( DWImagesType::Pointer &dwimage); 
  void GetComponentFromDWI( DWImagesType::Pointer dwimage, unsigned int n);
  void LoadROIS();
  void EstimateTensor( DWImagesType::Pointer dwimage, TensorImageType::Pointer &tensor, InputImageType::Pointer &t2image, bool mask);
  void ApplyMask( TensorImageType::Pointer &tensor);

  void verGlifos(int);
  void verGlifosTract();
//  void verGlifos(int,int,int,int,int,int,int);
  void borrarGlifos(int);
  void cambiarOpacidad(float);
  void glifosActivos();
  void imagenActiva(int,bool);
  void actualizarGlifosDTI();

  void verGlifosStrain();
  void verGlifosStrain(int,int);
  void borrarGlifosStrain();
  void ViewStrainSlice3D();
  void ViewStrainSlice3D(int);

  void SetScalarRangeStrain(double,double);
  void generarTensoresCilindro();
  void probarEsfuerzo();
  void probarEsfuerzoDeform();
  void tensorEsfuerzo2D();

  FloatImageType::IndexType     m_seedPosition;
  NodeType                      m_node;
  NodeContainer::Pointer        m_seeds;
  std::vector <std::vector<FloatImageType::IndexType> >  m_points_out;
  std::vector <std::vector<InputPointType > >  m_points_out_float;

  geodesicPath3DType::Pointer        m_geodesicPath3D;
  MyColorViewerType* ImageColorViewer;
  //Viewer3D* ImageViewer3D;
  UCharImageType* m_imageOverlay;
  Fl_Browser *m_modeldataBrowser;
  Fl_Browser *m_DWIdataBrowser;
  Fl_Browser *m_strainTensorDataBrowser;
  Fl_Value_Slider **sliceNumberSlider;
	
	typedef VolumesContainer<DataModelElementType>   VectorOfModelType;
	typedef VolumesContainer<DataElementType>        VectorOfDataType;
	typedef VolumesContainer<DataTensorElementType>  VectorOfTensorDataType;
	typedef VolumesContainer<DataDWIElementType>  VectorOfDWIDataType;
	typedef VolumesContainer<DataSTElementType>  VectorOfSTDataType;
	
	void SetVectorModelData(void* VectorModelData);
	void SetVectorTensorData(void* VectorTensorData);
	void SetVectorDWIData(void* VectorDWIData);
	void SetVectorData(void* VectorData);
	void SetVectorSTData (void* VectorSTData);
	
	void SetHasBeenRegistered(bool value){
		m_HasBeenRegistered=value;
	}

	void SetHasBeenFiberTracked(bool value){
		m_HasBeenFiberTracked=value;
	}
	
private: 
  TractographyType::Pointer			m_Tractography;
  DWIFilterType::Pointer            m_dwifilter;
  DTIEstimatorType::Pointer         m_dtiEstimator;
  FaTensorRegistrationType::Pointer m_register;
  
  VectorOfModelType*      m_VectorModelData;
  VectorOfDataType*       m_VectorData;
  VectorOfTensorDataType* m_VectorTensorData;
  VectorOfDWIDataType* m_VectorDWIData;
	Fl_Double_Window  *panelMessage;
	Fl_Text_Display *m_OrientationText,*m_OriginText; 
	Fl_Value_Slider *min_scalar_range,*max_scalar_range; 
  std::list<ClickPoint> clicked_points;

  VectorOfSTDataType*	m_VectorSTData;

  vtkActor		*m_activeActorX;
  vtkPolyData		*m_activeGlyphX;
  vtkActor		*m_activeActorY;
  vtkPolyData		*m_activeGlyphY;
  vtkActor		*m_activeActorZ;
  vtkPolyData		*m_activeGlyphZ;
  vtkScalarBarActor	*m_scalarBar;

  vtkPoints		*m_puntosTract;
  vtkActor		*m_tractActor;

  vtkActor		*m_activeActorStrain;
  vtkScalarBarActor	*m_scalarBarStrain;

  double rangeStrainMin;
  double rangeStrainMax;
  vtkFloatArray		*m_deformValues;


  vtkImageData 		*m_activeStrainImage;

  bool m_planoActivoX;
  bool m_planoActivoY;
  bool m_planoActivoZ;
  bool m_tractActiva;

  int tiempo;

  UCharImageType::Pointer	m_SeedRegionsImage;
  UCharImageType::Pointer	m_TractRegionsImage;
  ParametersType			m_Parameters;
  bool						m_HasBeenRegistered;
  bool						m_HasBeenFiberTracked;
  

};

#endif




