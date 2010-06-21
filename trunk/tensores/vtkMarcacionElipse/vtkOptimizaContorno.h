//////////////////////////////////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit		        				//
//	Modulo:		vtkOptimizaContorno.h	         					//
//	Descripcion:	Gestor de optimizacion	    	         				//
//	Fecha:		2005/04/04	     							//
//	Lenguaje:	C++						      			//
//	Autor:		Lucilio Cordero Grande	      						//
//			ETSI Telecomunicacion, Universidad de Valladolid			//
//			Campus Miguel Delibes, Camino del Cementerio, s/n			//
//			e-mail: lcorgra@lpi.tel.uva.es		      				//
//	Basado en un programa de:	Rocio Sanchez Fernandez	                		//
//					ETSI Telecomunicacion, Universidad de Valladolid	//
//					Campus Miguel Delibes, Camino del Cementerio, s/n	//
//					e-mail: rsanfer@atenea.lpi.tel.uva.es                   //
//////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef __vtkOptimizaContorno_h
#define __vtkOptimizaContorno_h

#include "vtkMarcacionElipsoideConfigure.h"
#include "vtkObject.h"
#include "vtkStructuredPoints.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include <math.h>

//class  VTK_MARCACIONELIPSOIDE_EXPORT vtkOptimizaContorno : public vtkObject
class  vtkOptimizaContorno : public vtkObject
{
 public:
   
   static vtkOptimizaContorno *New();
   const char *GetClassName() {return "vtkOptimizaContorno";};
   vtkTypeMacro(vtkOptimizaContorno, vtkObject);
   void PrintSelf(ostream& os, vtkIndent indent);
   
   /*Ns=barridos de SA*/
   vtkSetMacro(Ns,int);
   vtkGetMacro(Ns,int);
   
   vtkSetStringMacro(Ruta);
   vtkGetStringMacro(Ruta);
   
   vtkSetVector3Macro(DimMalla,int);
   vtkGetVector3Macro(DimMalla,int);
   
   vtkSetVector3Macro(DimZ,int);
   vtkGetVector3Macro(DimZ,int);
   
   vtkSetVector3Macro(PeriodMalla,int);
   vtkGetVector3Macro(PeriodMalla,int);
   
   vtkSetMacro(NumEntidades,int);
   vtkGetMacro(NumEntidades,int);
      
   vtkSetMacro(DimensionalidadEstados,int);
   vtkGetMacro(DimensionalidadEstados,int);
   
   vtkSetMacro(DimensionalidadMalla,int);
   vtkGetMacro(DimensionalidadMalla,int);
   
   vtkSetMacro(Independencia,int);
   vtkGetMacro(Independencia,int);
   
   vtkSetMacro(OrdenMalla,float);
   vtkGetMacro(OrdenMalla,float);
   
   vtkSetMacro(OrdenEstados,float);
   vtkGetMacro(OrdenEstados,float);
   
   vtkSetMacro(V,int);
   vtkGetMacro(V,int);
   
   vtkSetMacro(K,int);
   vtkGetMacro(K,int);
   vtkSetMacro(Q,int);
   vtkGetMacro(Q,int);
   
   vtkSetMacro(Super,int);
   vtkGetMacro(Super,int);
   
   vtkSetMacro(UsimagTool,int);
   vtkGetMacro(UsimagTool,int);

   vtkStructuredPoints *DevuelveRho();

   
   void EstableceLE(vtkStructuredPoints *,int,int);
   void EstableceRhoMed(vtkStructuredPoints *);
   void EstableceRho(vtkStructuredPoints *);
   void Llamada(int);
   void ConstruyeVecindario();
   
 protected:
   vtkOptimizaContorno();
   ~vtkOptimizaContorno();
   
   vtkPoints *Puntos;
   vtkPoints *PuntObservacion;
   vtkPoints *PuntParam;

   vtkPoints *Vecindario;
   vtkIntArray *NumClNivel;
   vtkCellArray *Cliques;
   vtkCellArray *OrdenCliques;
   vtkCellArray *CatCliques;
   vtkStructuredPoints *CampoMedio[2];
   vtkStructuredPoints *LE[2][20];
   vtkStructuredPoints *EnergiaNormalizada[20];
   vtkStructuredPoints *Potencia[20];
   
   vtkStructuredPoints *Rho;
   vtkStructuredPoints *RhoBordes;
   vtkStructuredPoints *RhoMarMed;
   vtkStructuredPoints *Zonas;
   vtkStructuredPoints *Parametros;
   vtkStructuredPoints *MetricaCliques;

   vtkFloatArray *Multiplicadores;
   vtkDoubleArray *EnE[2];
   vtkDoubleArray *Pg;
   //Alfa, en el caso de un campo no homogeneo (inclusion de la coordenada 
   //en profundidad) deberia tener mas componentes).
   vtkDoubleArray *AlfaPot;
   vtkDoubleArray *BetaPot;
   vtkFloatArray *Deformaciones;

   
//   void Subredecilla(int);
//   float Energia(int,int,int,bool,float *);
   void Subredecilla(int,int,int);
   void ActualizaRhoBordes(int,int,int,int,int);
   float Energia(int *,int,int,int,vtkFloatArray *);
   float Potencial(float *,int);

   void ConvierteRhoBordes();
   void ConvierteRho();
   void Ejecuta();
   void OptimizaParametros();
   void CalculaEnergiaNormalizada();
   void Monitoriza();
   void ConstruyeMetricaCliques();
   
   /*barridos de SA*/
   int Ns;
   
   /*barridos de Metropolis-Hastings*/
   int Nsparam;
   
   /*barridos SA entre cada barrido Metropolis-Hastings*/
   int Nsinter;
   
   int DimMalla[3];             //Dimensiones del campo de deformaciones
   int DimZ[3];
   int PeriodMalla[3];          //Periodicidad de cada una de las dimensiones (1=Ciclico,0=NoCiclico)
   int NumEntidades;            //Numero de campos de deformaciones a considerar.
   int DimensionalidadEstados;  //Numero de dimensiones del espacio de posibles deformaciones.
   int DimensionalidadMalla;    //Numero de dimensiones no singulares de la malla.
   int Independencia;           //(1=Considera solo cliques unidimensionales/0=Considera cliques bidimensionales).
   //int OrdenMalla;              //Orden del vecindario.
   //int OrdenEstados;            //Orden del vecindario de deformacion.
   float OrdenMalla;
   float OrdenEstados;
   bool MFA;                    //Indica si se va a hacer el Mean Field Annealing.
   
   /*cardinalidad del espacio de estados*/
   int K;
   int escribe;
   float Beta;
   float FactConv;
   int P;
   int Q;
   int O;
   int Super;
   /*número de parámetros por cada contorno*/
   int L;
   /*número de términos de verosimilitud*/
   int V;
   bool Leido;
   char *Ruta;
   int Ord;
   double Suma; //Modulo del vector de parametros.
   int UsimagTool;
};

#endif
