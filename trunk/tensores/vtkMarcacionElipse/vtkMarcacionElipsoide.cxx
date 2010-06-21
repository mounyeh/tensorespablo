/////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit
//	Módulo:	vtkMarcacionElipsoide.cxx
//	Descripción:	Genera un elipsoide para interaccion con 3D slicer
//	Fecha:	2004/08/31
//	Última modificación: 2006/05/22
//	Lenguaje:	C++
//	Autor:	Lucilio Cordero Grande
//		ETSI Telecomunicacion, Universidad de Valladolid
//		Campus Miguel Delibes, Camino del Cementerio, s/n
//		e-mail: lcorgra@atenea.lpi.tel.uva.es
//
//////////////////////////////////////////////////////////////////////

#include "vtkMarcacionElipsoide.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPointData.h"
#include "vtkLookupTable.h"

vtkMarcacionElipsoide* vtkMarcacionElipsoide::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMarcacionElipsoide");
  if(ret)
    {
    return (vtkMarcacionElipsoide*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMarcacionElipsoide;
}

vtkMarcacionElipsoide::vtkMarcacionElipsoide()
{
   this->NomFichMar=NULL;
   this->SetNomFichMar("PuntosMarcados.vtk");
   this->BarraColor=vtkScalarBarActor::New();
}

vtkMarcacionElipsoide::~vtkMarcacionElipsoide()
{
}

//Genera un elipsoide a partir de sus ejes principales.
vtkPolyData *vtkMarcacionElipsoide::GeneraElipsoide(vtkPoints *Ejes)
{
   vtkPolyData *Poli=vtkPolyData::New();

   Poli->SetPoints(Ejes);
   Poli->Update();

   return Poli;
}

void vtkMarcacionElipsoide::GuardarMarcacionInicial(vtkPoints *Ej)
{
   int j;
   char nombre[200];
   
   vtkPolyDataWriter *Escritor=vtkPolyDataWriter::New();
   vtkPolyData *Datos=vtkPolyData::New();

   j=sprintf(nombre,"%s.Punt.vtk",this->NomFichMar);
   //j=sprintf(nombre,"%s",this->NomFichMar);
   Datos->SetPoints(Ej);
   Escritor->SetFileName(nombre);
   Escritor->SetFileTypeToASCII();
   Escritor->SetInput(Datos);
   Escritor->Write();
}

vtkPoints *vtkMarcacionElipsoide::CargarMarcacionInicial()
{
   int j;
   char nombre[200];
   
   vtkPoints *Ej=vtkPoints::New();
   vtkPolyDataReader *Lector=vtkPolyDataReader::New();
   vtkPolyData *Datos=vtkPolyData::New();

   j=sprintf(nombre,"%s.Punt.vtk",this->NomFichMar);
   //j=sprintf(nombre,"%s",this->NomFichMar);
   Lector->SetFileName(nombre);
   Datos=Lector->GetOutput();
   Datos->Update();
   Ej=Datos->GetPoints();

   return Ej;
}

vtkActor *vtkMarcacionElipsoide::DibujaModelo(vtkPolyData *poli,int tipo,int col)
{
   int i;
   
   vtkPolyDataMapper *map_elip=vtkPolyDataMapper::New();
   vtkActor *actor_elip=vtkActor::New();
   vtkProperty *color=vtkProperty::New();
   vtkLookupTable *mapeoescalar=vtkLookupTable::New();

   mapeoescalar->SetTableRange(this->Min[col],this->Max[col]);
   if (col==0)
     {
	color->SetColor(1,0,0);
	i=poli->GetPointData()->SetActiveScalars("CurvGauss");
	for (i=0;i<256;i++)
	  mapeoescalar->SetTableValue(i,i/255.0,.25,.75,1);
	this->BarraColor->SetTitle("Curvat. gaussi.");
     }
   if (col==1)
     {
	color->SetColor(0,1,0);
	i=poli->GetPointData()->SetActiveScalars("CurvMedia");
	for (i=0;i<256;i++)
	  mapeoescalar->SetTableValue(i,i/255.0,.75,.25,1);
	this->BarraColor->SetTitle("Curvatura media");
     }
   if (col==2)
     {
	color->SetColor(0,0,1);
	i=poli->GetPointData()->SetActiveScalars("IndForma");
	for (i=0;i<256;i++)
	  mapeoescalar->SetTableValue(i,i/255.0,0,1,1);
	this->BarraColor->SetTitle("Indice de forma");
     }
   if (col==3)
     {
	color->SetColor(1,1,0);
	i=poli->GetPointData()->SetActiveScalars("IndCurva");
	for (i=0;i<256;i++)
	  mapeoescalar->SetTableValue(i,i/255.0,1,0,1);
	this->BarraColor->SetTitle("Indice de curv.");
     }
   if (col==4)
     color->SetColor(1,0,1);
   if (col==5)
     color->SetColor(0,1,1);
   
   map_elip->SetLookupTable(mapeoescalar);
   map_elip->SetScalarRange(this->Min[col],this->Max[col]);
   this->BarraColor->SetLabelFormat("%.2f             ");
   this->BarraColor->SetLookupTable(mapeoescalar);
   map_elip->SetColorModeToMapScalars();
   if (tipo==1)
     map_elip->ScalarVisibilityOn();
   map_elip->SetInput(poli);
   map_elip->GlobalImmediateModeRenderingOn();
   map_elip->Update();
   actor_elip->SetMapper(map_elip);
   color->SetOpacity(0.5);
   color->SetDiffuse(0.0);
   color->SetSpecular(0.0);
   color->SetAmbient(0.9);
   actor_elip->SetProperty(color);
   
		     mapeoescalar->Delete();
   map_elip->Delete();
   color->Delete();
   return actor_elip;
}

vtkScalarBarActor *vtkMarcacionElipsoide::DevuelveBarra()
{
   this->BarraColor->SetNumberOfLabels(5);
   
   return this->BarraColor;
}

void vtkMarcacionElipsoide::PrintSelf(ostream& os, vtkIndent indent)
{
}
