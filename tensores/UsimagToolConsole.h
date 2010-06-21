/*=========================================================================

  Program:   UsimagToolConsole
  Language:  C++
  Date:      5-07-2007
  Version:   1.0

  Copyright (c) 2007 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/

#include "UsimagToolGUI.h"
#include <itkObject.h>
#include <string>
#include <itkImageToImageFilter.h>

namespace itk {
  template<class TObject> class MemberCommand;
  class EventObject;
}

class UsimagToolConsole : public UsimagToolGUI
{
public:
  typedef itk::MemberCommand<UsimagToolConsole> CommandType;
  typedef itk::SmartPointer<CommandType> CommandPointer;
  typedef itk::ImageToImageFilter<InputImageType,InputImageType> itkFilter;
  
  UsimagToolConsole();
  virtual ~UsimagToolConsole();
  void Save(void);
  void SaveModel(void);
  void SaveTensor(void);
  void SaveTensorVTK(void);
  void Save(const char*);
  void SaveOverlay(void);
  void SaveDWI(void);
  void Load(void);
  void LoadTensor(void);
  void LoadME(void);
  void Abrir(int filetype);
  void AbrirTensor(int filetype);
  void Quit(void);
  void LoadRaw(void);
  void LoadRaw(const char* filename);
  void LoadMeta(void); 
  void LoadNifti( void );
  void LoadVOL(void);
  void LoadColor(void);
  void LoadJPEG(void);
  void LoadVTK(void);
  void LoadTensorVTK(void);
  void LoadTensorNrrd(void);
  void LoadTIFF(void);
  void LoadPNG(void);
  void LoadDicom(void);
  void LoadKNN(void);
  void LoadGenericDWI(void);
  void LoadGenericDWIandProcess(void);
  void LoadDicomDWI(void);
  void LoadDicomDWIandProcess(void);
  void LoadNrrdDWI(void);
  void LoadNrrdDWIandProcess(void);
  void LoadDWI(void);
  void LoadDicomSeries(void);
  void LoadNrrd();
  void LoadModel(void);
  void LoadFibers(void);
  void LoadStrainTensor();
  void LoadStrainTensorNrrd();
  void WriteFibers(vtkPolyData* streamlines);
  void Distribuir(InputImageType::Pointer aux, int window, const char* filename);
  void DistribuirColor(ImageRGBType::Pointer aux, int window, const char* filename);
  void DistribuirTensor(TensorImageType::Pointer aux, int window, const char* filename);
  int check_orientation(InputImageType::DirectionType Dir);
  void OrientarImagen(InputImageType::Pointer &image, int orientacion);
  
  // Metodos 
  void ShowImageSrc(int window);
  void ShowImageColor(int window);
  void ViewSliceIn3DWindow(InputImageType::Pointer image, int orientation);
  void DoModel(int value1, int value2, int smooth, InputImageType::Pointer image);
  void OnProgress(itk::Object *source, const itk::EventObject &event);
  void OnFileInputChange(int mode);
  void OnFileFormatChange(int mode);
  void OnTipoProcesadoChange();
  int  DetermineFileFormatFromFileName(const char *testFile,int mode);
  void OnPixelTypeChange();
  void OnByteOrderChange();
  void OnSliceChange(unsigned int value, unsigned int window);
  void OnOrientationChange(unsigned int window);
  void ChangeViewerColorMode();
  void OnCheckButtonChange(bool t, int window);
  void OnImageModeChange();
  void ViewData(int, int);
  void LoadOverlay(int);	
  void InitializeOverlay(int window);
  void SetOverlay(int window);
  void ClearOverlay();
  void DeleteScalar( int );
  void DeleteTensor( int );
  void DeleteModel( int );
  void RenameScalar( int, const char* );
  void RenameTensor( int, const char* );
  void RenameModel( int, const char* );
  void GenerateNewData(InputImageType::Pointer aux);
  void CopyData(InputImageType::Pointer aux);
  void OnAmplify(Fl_Group* viewSplit,int window);
  void OnMinimize(int window,int x, int y);
  void UpdateIntensityW(int window, float val);
  void UpdateIntensityL(int window, float val);
  void directorio(const char* file,char* dir);
  void str_nombre(const char* file,char* nombre);
  void ViewMode3D();
  void ViewMode4x2D();
  void ViewMode3_1();
  void InfoImagen( void );
  

  void GenericFilter(itkFilter*);
  // Filtros 
  // Registration 
  void RegisterMLDFilter(void);
  //typedef itk::Object CommandType;
  //typedef itk::Object::Pointer CommandPointer;
  
  typedef InputImageType::PointType						InputPointType;
  typedef std::vector<InputPointType>					InputPointVectorType;
  std::vector<InputPointVectorType>						m_points_out_float;
  
  int                                                   m_filetype;

  //=========================================================
  /** Zooming the image in */
  void ZoomIn();
  void ZoomOut();
  //=========================================================
  
  
private:
  static const int m_FORMAT_COUNT = 10; 
  std::string m_FileFormatPattern[m_FORMAT_COUNT];
  static const int m_TENSOR_FORMAT_COUNT = 2; 
  std::string m_FileTensorFormatPattern[m_FORMAT_COUNT];
  int m_orientation;
  
};

