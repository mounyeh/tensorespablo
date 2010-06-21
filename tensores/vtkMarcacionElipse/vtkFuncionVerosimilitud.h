/////////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit								   //
//	Módulo:	vtkFuncionVerosimilitud.h								   //
//	Descripción:	Recoge métodos para estimar los valores de las dos //
//			funciones de verosimilitud empleadas en el algoritmo	   //
//	Fecha:	2004/09/09												   //
//	Lenguaje:	C++													   //
//	Autor:	Lucilio Cordero Grande									   //
//		ETSI Telecomunicacion, Universidad de Valladolid			   //
//		Campus Miguel Delibes, Camino del Cementerio, s/n			   //
//		e-mail: lcorgra@atenea.lpi.tel.uva.es						   //
/////////////////////////////////////////////////////////////////////////

#ifndef __vtkFuncionVerosimilitud_h
#define __vtkFuncionVerosimilitud_h

#include <vtkMarcacionElipsoideConfigure.h>
#include "vtkObject.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkStructuredPoints.h"
#include "vtkTIFFWriter.h"
#include "vtkMatrix4x4.h"
#include "vtkIntArray.h"

//class VTK_MARCACIONELIPSOIDE_EXPORT vtkFuncionVerosimilitud : public vtkObject
class vtkFuncionVerosimilitud : public vtkObject
{
public:
  static vtkFuncionVerosimilitud *New();
  vtkTypeMacro(vtkFuncionVerosimilitud, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  const char *GetClassName() {return "vtkFuncionVerosimilitud";};

  void IntroduceParam(vtkStructuredPoints *);
  void IntroduceCentro(vtkFloatArray *);
  void IntroduceImagen(vtkImageData *,int);
  void IntroduceMatriza(vtkMatrix4x4 *);
  void IntroducePlano(vtkMatrix4x4 *);
    
  /*limites*/
  vtkSetVector6Macro(limite,double);
  vtkGetVector6Macro(limite,double);
  
  vtkSetMacro(Ng,int);
  vtkGetMacro(Ng,int);

  vtkSetMacro(J,int);
  vtkGetMacro(J,int);
   
   vtkSetMacro(UsimagTool,int);
   vtkGetMacro(UsimagTool,int);
  
  vtkSetMacro(P,int);
  vtkGetMacro(P,int);
  
  vtkSetMacro(K,int);
  vtkGetMacro(K,int);
  
  //vtkSetMacro(M,int);
  vtkGetMacro(Mbis,int);
  
  //vtkSetMacro(N,int);
  vtkGetMacro(Nbis,int);
  
  vtkGetMacro(Z,int);
 
  vtkSetMacro(drmax,float);
  vtkGetMacro(drmax,float);
   
   vtkSetStringMacro(Ruta);
   vtkGetStringMacro(Ruta);
 
   void Ejecuta();
   int CalculaBetaGamma(int);
   void CalculaIf();
   void CalculaG();
   void CalculaLB(int,int);
   void CalculaLI(int,int);
   //Metodos de acceso externo
   void EscribeImagen();
   vtkStructuredPoints *DevuelveLR();
   vtkIntArray *DevuelveMascaraElipsoide();

protected:

  vtkFuncionVerosimilitud();
  ~vtkFuncionVerosimilitud();

   float mu;
   float g;
   float drmax;
   int Ng;
   
   int salto;
   int J;
   int P;
   int K;
   int M;
   int N;
   int Z;
   int Mbis;
   int Nbis;
   int UsimagTool;
   
   double limite[6];
   int limites[6];
   float *dr;
   int *modmed;
   char *Ruta;

  vtkStructuredPoints *Beta;
  vtkStructuredPoints *Gamma;
  vtkStructuredPoints *Mask;
  vtkMath *Instrumento;
  vtkStructuredPoints *mod;  //La primera componente se aplica a intensidad(gamma) y la segunda a gradiente(beta).
  vtkStructuredPoints *LR;  //Verosimilitud de intensidad en componente 0 y de borde en componente 1.
  vtkImageData *ImEnt[2];  //La 0 es la normal y la 1 la filtrada.
  vtkImageData *Sub;
  vtkImageData *If;
  vtkImageData *G;
  vtkStructuredPoints *Param;
  vtkFloatArray *Centro;
  vtkIntArray *MascaraElipsoide;

  vtkMatrix4x4 *MatrizRASaIJK;
  vtkMatrix4x4 *MatrizIJKaRAS;
  vtkMatrix4x4 *Plano;

};

#endif
