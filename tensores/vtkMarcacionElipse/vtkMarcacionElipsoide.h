/////////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit								   //
//	Módulo:		vtkMarcacionElipsoide.h								   //
//	Descripción:	Genera un elipsoide para interaccion con 3D slicer //
//	Fecha:		2004/08/31											   //
//	Lenguaje:	C++													   //
//	Autor:		Lucilio Cordero Grande								   //
//			ETSI Telecomunicacion, Universidad de Valladolid		   //
//			Campus Miguel Delibes, Camino del Cementerio, s/n		   //
//			e-mail: lcorgra@atenea.lpi.tel.uva.es					   //
/////////////////////////////////////////////////////////////////////////



// .NOMBRE vtkMarcacionElipsoide - genera un elipsoide que supone el primer ajuste
// .SECCION Descripcion
// A partir de 6 puntos que suponen los extremos de los ejes principales genera un elipsoide

#ifndef __vtkMarcacionElipsoide_h
#define __vtkMarcacionElipsoide_h

#include "vtkMarcacionElipsoideConfigure.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"


//class VTK_MARCACIONELIPSOIDE_EXPORT vtkMarcacionElipsoide : public vtkObject
class vtkMarcacionElipsoide : public vtkObject
{
public:
  static vtkMarcacionElipsoide *New();
  vtkTypeMacro(vtkMarcacionElipsoide, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

   vtkSetStringMacro(NomFichMar);
   vtkGetStringMacro(NomFichMar);
   
   vtkSetVectorMacro(Max,double,4);
   vtkGetVectorMacro(Max,double,4);
   
   vtkSetVectorMacro(Min,double,4);
   vtkGetVectorMacro(Min,double,4);
  
  vtkPolyData *GeneraElipsoide(vtkPoints *);
  vtkActor *DibujaModelo(vtkPolyData *,int,int);
  void GuardarMarcacionInicial(vtkPoints *);
  vtkPoints *CargarMarcacionInicial();
  vtkScalarBarActor *DevuelveBarra();

protected:
  vtkMarcacionElipsoide();
  ~vtkMarcacionElipsoide();
   
   char *NomFichMar;
   vtkScalarBarActor *BarraColor;
   double Max[4];
   double Min[4];
};
#endif
