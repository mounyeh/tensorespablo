/*=========================================================================

  Program:   UsimagToolConsoleSegmentation
  Language:  C++
  Date:      5-07-2007
  Version:   1.0

  Copyright (c) 2007 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#include "UsimagToolConsoleSegmentation.h"
#include "vtkITKUtility.h"
#include "vtkPolyDataWriter.h"


UsimagToolConsoleSegmentation::UsimagToolConsoleSegmentation() {
  m_GradientFilter =  GradientFilterType::New();
}

UsimagToolConsoleSegmentation::~UsimagToolConsoleSegmentation() {
}


void UsimagToolConsoleSegmentation::GradientFilter( unsigned int threads, InputImageType* inputImage, InputImageType* outputImage )
{

  std::cout << "Haciendo Gradiente " << std::endl;  
  //CommandPointer callback = CommandType::New();  
  //callback->SetCallbackFunction(this,&UsimagToolConsoleSegmentation::OnProgress);

  m_GradientFilter->SetNumberOfThreads(threads);
  m_GradientFilter->SetInput( inputImage );
  //m_GradientFilter->AddObserver( itk::ProgressEvent(), callback);
  m_GradientFilter->Update();

  if (m_destino->value() == 0) {
    outputImage = m_GradientFilter->GetOutput();
  } else {
    inputImage[m_destino->value()-1] = m_GradientFilter->GetOutput();
    this->ShowImageSrc(m_destino->value()-1);
  }

  // Ponemos la barra de progreso a cero  
  progressSlider->value(0);
  m_ProgressCounter->value(0);
}

void UsimagToolConsoleSegmentation::MEFileInputChange() {
  const char *text = ME_FilePageBrowser->value();
  if (text != NULL && strlen(text) > 0) {    
    ME_filename = text;
  }
}

void UsimagToolConsoleSegmentation::MarcacionElipseInit() {
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

void UsimagToolConsoleSegmentation::CalculoModeloVerosimilitud() {

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
  m_VTKexporter->SetInput(m_InputImage[m_activeinput]);
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
     tran->SetElement(i,3,(int)m_InputImage[m_activeinput]->GetLargestPossibleRegion().GetSize()[i]/2);
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

void UsimagToolConsoleSegmentation::CalculoOptimizacionSA() {
   
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

