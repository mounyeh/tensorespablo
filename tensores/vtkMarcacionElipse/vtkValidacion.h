/////////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit								   //
//	Módulo:	vtkValidacion.h								   //
//	Descripción:	Recoge métodos para obtener los cortes con cada imagen del modelo
//			obtenido, con objetivo de compararlos con los resultados médicos
//	Fecha:	2005/09/09												   //
//	Lenguaje:	C++													   //
//	Autor:	Lucilio Cordero Grande									   //
//		ETSI Telecomunicacion, Universidad de Valladolid			   //
//		Campus Miguel Delibes, Camino del Cementerio, s/n			   //
//		e-mail: lcorgra@troya.lpi.tel.uva.es						   //
/////////////////////////////////////////////////////////////////////////

#ifndef __vtkValidacion_h
#define __vtkValidacion_h

#include "vtkMarcacionElipsoideConfigure.h"
#include "vtkObject.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkMatrix4x4.h"
#include "vtkPoints.h"

//class VTK_MARCACIONELIPSOIDE_EXPORT vtkValidacion : public vtkObject
class vtkValidacion : public vtkObject
{
public:
  static vtkValidacion *New();
  vtkTypeMacro(vtkValidacion, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  const char *GetClassName() {return "vtkValidacion";};

  //Metodos de acceso externo

void IntroduceMatriza(vtkMatrix4x4 *);
void IntroducePlantilla(vtkPolyData *);
void GeneraValores();
vtkPolyData *Ejecuta();

protected:

  vtkValidacion();
  ~vtkValidacion();



   double Origen[3];
   double Incremento;

  vtkFloatArray *OrigenImag;
  vtkFloatArray *UnitXImag;
  vtkFloatArray *UnitYImag;
  vtkPolyData *Plantilla;
  vtkPoints *ConjPunt;
  
  vtkMatrix4x4 *MatrizRASaIJK;
  vtkMatrix4x4 *MatrizIJKaRAS;

};

#endif
