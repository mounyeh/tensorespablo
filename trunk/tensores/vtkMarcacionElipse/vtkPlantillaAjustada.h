//////////////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit										//
//	Módulo:	vtkPlantillaAjustada.h											//
//	Descripción:	A partir de una superficie que aproxima el contorno 3D  //
//				genera una serie de parámetros útiles en nuestro algoritmo. //
//	Fecha:	2004/09/04														//
//	Lenguaje:	C++															//
//	Autor:	Lucilio Cordero Grande											//
//		ETSI Telecomunicacion, Universidad de Valladolid					//
//		Campus Miguel Delibes, Camino del Cementerio, s/n					//
//		e-mail: lcorgra@atenea.lpi.tel.uva.es								//
//////////////////////////////////////////////////////////////////////////////


// .NOMBRE vtkPlantillaAjustada - genera centro, rayos y radios para cada slice.
// .SECCION Descripcion
// Lleva a cabo la parametrización que se utiliza para la superficie externa de órgano 
// en nuestro algoritmo.

#ifndef __vtkPlantillaAjustada_h
#define __vtkPlantillaAjustada_h


#include "vtkMarcacionElipsoideConfigure.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkPolyData.h"
#include "vtkStructuredPoints.h"
#include "vtkIntArray.h"

//class VTK_MARCACIONELIPSOIDE_EXPORT vtkPlantillaAjustada : public vtkObject
class vtkPlantillaAjustada : public vtkObject
{
 public:
   static vtkPlantillaAjustada *New();
   vtkTypeMacro(vtkPlantillaAjustada, vtkObject);
   void PrintSelf(ostream& os, vtkIndent indent);
   //void GeneraParametros (vtkPlane *);
   //void IntroducePlantilla (vtkPolyData *);
   void ObtieneParametros(vtkPolyData *);
   
   vtkFloatArray *ObtieneCentro();
   vtkStructuredPoints *ObtieneParam();
  
   void IntroducePlano(vtkMatrix4x4 *);
   void EstableceMascaraElipsoide(vtkIntArray *);
   vtkPolyData *ConstruyeModelo(vtkStructuredPoints *);
   vtkStructuredPoints *GeneraRhoNulo();
   vtkPolyData *LeeModelo();
  
   /*J=número de rayos*/
   vtkSetMacro(J,int);
   vtkGetMacro(J,int);

   /*K=número de puntos de la deformación*/
   vtkSetMacro(K,int);
   vtkGetMacro(K,int);
  
   /*Slth=distancia entre slices*/
   vtkSetMacro(Slth,float);
   vtkGetMacro(Slth,float);
   
   vtkSetMacro(UsimagTool,int);
   vtkGetMacro(UsimagTool,int);
  
   /*P=número de slices de interés*/
   vtkGetMacro(P,int);
  
   /*Máxima desviación radial*/
   vtkSetMacro(drmax,float);
   vtkGetMacro(drmax,float);
  
   /*Limites de la caja de inclusión del ajuste inicial*/
   vtkSetVector6Macro(limite,double);
   vtkGetVector6Macro(limite,double);
   
   vtkSetVectorMacro(Max,double,4);
   vtkGetVectorMacro(Max,double,4);
   
   vtkSetVectorMacro(Min,double,4);
   vtkGetVectorMacro(Min,double,4);
   
   vtkSetStringMacro(FichModelo);
   vtkGetStringMacro(FichModelo);

protected:
   vtkPlantillaAjustada();
   ~vtkPlantillaAjustada();
   vtkFloatArray *Centro;
   vtkStructuredPoints *Param; //En la componente 0 esta el angulo y en la 1 el radio.
   vtkMath *Instrumento;
   vtkMatrix4x4 *Plano;
   vtkIntArray *MascaraElipsoide;
   vtkPolyData *Rinon;

   int K;
   int J;
   int P;
   int salto;
   int UsimagTool;
   float drmax;
   float Slth;
   double limite[6];
   double Max[4];
   double Min[4];
   char *FichModelo;
};
#endif
