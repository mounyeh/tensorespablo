////////////////////////////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit						  //
//	Modulo:	vtkOptimizaContorno.cxx			         			  //
//	Descripcion:	Gestor de optimizacion de la posicion del contorno                //
//			mediante el algoritmo SA.				      	  //
//	Fecha:	2006/01/10					     			  //
//	Lenguaje:	C++		     						  //
//	Autor:	Lucilio Cordero Grande							  //
//		ETSI Telecomunicacion, Universidad de Valladolid			  //
//		Campus Miguel Delibes, Camino del Cementerio, s/n			  //
//		e-mail: lcorgra@lpi.tel.uva.es			       	              	  //
//	Basado en un programa de:	Rocio Sanchez Fernandez				  //
//					ETSI Telecomunicacion, Universidad de Valladolid  //
//					Campus Miguel Delibes, Camino del Cementerio, s/n //
//					e-mail: rsanfer@atenea.lpi.tel.uva.es		  //
////////////////////////////////////////////////////////////////////////////////////////////


//En la tesis de Marcos aparece una ecuacion que induce a pensar que para 
//la independencia entre puntos de la malla habria que tomar 5 unidades en 
//vez de 3 unidades.

/*
 * Para el trabajo de investigacion a llevar a cabo para justificar el segundo a√±o, resultaria interesante hacer esta clase cutre una clase en condiciones, permitiendo fijar entre otras cosas:
 * El tama√±o del vecindario.
 * Los cliques no nulos (y las anisotropias).
 * Los terminos de verosimilitud de entrada.
 * El structured points de salida.
 * A la entrada tenemos un vtkStructuredPoints de n componentes y a la salida un vtkStructuredPoints de una sola componente.
 * Por tanto habria que hacerle heredar de vtkStructuredPointsToStructuredPointsFilter.
 * Otros parametros serian el caracter ciclico o no del vecindario, la estimacion de parametros supervisada o no (e incluso semisupervisada) y el numero de campos a estimar.
 * Este dise√±o general habria que hacerlo lo antes posible.
 */


#include "vtkOptimizaContorno.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkStructuredPointsWriter.h"


vtkOptimizaContorno* vtkOptimizaContorno::New()
{
   // First try to create the object from the vtkObjectFactory
   vtkObject* ret = vtkObjectFactory::CreateInstance("vtkOptimizaContorno");
   if(ret)
     {
	return (vtkOptimizaContorno*)ret;
     }
   // If the factory was unable to create the object, then create it here.
   return new vtkOptimizaContorno;
}

vtkOptimizaContorno::vtkOptimizaContorno()
{
   int i;
   /*valores para SA*/
   this->Q=0;
   this->P=0;
   this->O=0;
   this->MFA=0;
   this->Ns=150;
   //Lo reduzco a una dÈcima parte para desprenderme de carga computacional:
   //this->Nsparam=1000;	//Quiz· habrÌa que bajarlo; est· asÌ para la supervisada (50 en la no supervisada).
   //	this->Nsparam=10;
   this->Nsparam=10;    //Reduciendo el tiempo de procesado para primeras pruebas
   this->Nsinter=25;
   this->FactConv=24.0;
   this->Suma=1.0;

   this->escribe=1;
   this->Leido=0;
   this->UsimagTool=0;
   
   this->Rho=vtkStructuredPoints::New();
   this->RhoBordes=vtkStructuredPoints::New();
   this->Zonas=vtkStructuredPoints::New();
   this->MetricaCliques=vtkStructuredPoints::New();
   
   for (i=0;i<2;i++)
     this->EnE[i]=vtkDoubleArray::New();
   this->Pg=vtkDoubleArray::New();
   this->Parametros=vtkStructuredPoints::New();
   this->Puntos=vtkPoints::New();
   this->PuntObservacion=vtkPoints::New();
   this->PuntParam=vtkPoints::New();
   this->Vecindario=vtkPoints::New();
   this->Cliques=vtkCellArray::New();
   this->OrdenCliques=vtkCellArray::New();
   this->CatCliques=vtkCellArray::New();
   this->NumClNivel=vtkIntArray::New();
   this->AlfaPot=vtkDoubleArray::New();
   this->BetaPot=vtkDoubleArray::New();
   this->Multiplicadores=vtkFloatArray::New();
   this->Deformaciones=vtkFloatArray::New();

   this->Puntos->SetDataTypeToFloat();
   this->PuntObservacion->SetDataTypeToFloat();
   this->PuntParam->SetDataTypeToFloat();
   
   this->Ruta=NULL;
//   Caso renal:
   this->SetRuta("");
   //Caso cardiaco:
//   this->SetRuta("Modules/vtkSegmentaMiocardio/datos/mra/se1");
}

vtkOptimizaContorno::~vtkOptimizaContorno()
{
   int i,j;
   
   this->Rho->Delete();
   this->RhoBordes->Delete();
   this->Zonas->Delete();
   this->Pg->Delete();
   this->Parametros->Delete();
   this->Puntos->Delete();
   this->PuntObservacion->Delete();
   this->PuntParam->Delete();
   
   for (i=0;i<this->NumEntidades;i++)
     this->EnE[i]->Delete();
     {
	for (j=0;j<this->DimMalla[2];j++)
	  this->LE[i][j]->Delete();
     }
   for (j=0;j<this->DimMalla[2];j++)
     {
	this->Potencia[j]->Delete();
	this->EnergiaNormalizada[j]->Delete();
     }
   this->RhoMarMed->Delete();
}	

void vtkOptimizaContorno::Llamada(int manera)
{		
   int k,kbis,j,kmin,n,naux,z,a,t,u,inda[3],indb[3];
   float min,aux,*tuplas;
   double punto[3],paraux[22];
   char nombre[200];
   char nomb[200];

   vtkStructuredPointsWriter *Escritora=vtkStructuredPointsWriter::New();
   vtkPolyDataWriter *Escritorb=vtkPolyDataWriter::New();
   //vtkPolyDataWriter *Escritorc=vtkPolyDataWriter::New();
   vtkPolyDataReader *Lector=vtkPolyDataReader::New();
   vtkPolyData *Auxiliar=vtkPolyData::New();
   vtkPolyData *Poli=vtkPolyData::New();
   vtkStructuredPoints *RhoEnfr=vtkStructuredPoints::New();
   vtkStructuredPointsWriter *EscribeRho=vtkStructuredPointsWriter::New();
   
   tuplas=new float[this->DimensionalidadEstados];
   
   this->Deformaciones->Initialize();
   
   this->Deformaciones->SetNumberOfComponents(this->DimensionalidadEstados);
   this->Deformaciones->SetNumberOfTuples(this->K);
   
   kbis=(int)floor(sqrt(this->OrdenEstados));
   
   if (this->DimensionalidadEstados==1)
     {
	for (k=0;k<this->K;k++)
	  this->Deformaciones->SetValue(k,-(this->K-1)/2+k);
     }
   else
     {
	j=0;
	for (k=0;k<3;k++)
	  inda[k]=0;
	for (k=0;k<this->DimensionalidadEstados;k++)
	  {
	     inda[k]=1;
	     tuplas[k]=0;
	  }
	this->Deformaciones->SetTuple(0,tuplas);
	j++;
	for (indb[2]=-kbis*inda[2];indb[2]<=kbis*inda[2];indb[2]++)
	  {
	     for (indb[1]=-kbis*inda[1];indb[1]<=kbis*inda[1];indb[1]++)
	       {
		  for (indb[0]=-kbis*inda[0];indb[0]<=kbis*inda[0];indb[0]++)
		    {
		       aux=pow((float)indb[0],2)+pow((float)indb[1],2)+pow((float)indb[2],2);
		       if (aux!=0 && aux<=this->OrdenEstados)
			 {
			    for (k=0;k<this->DimensionalidadEstados;k++)
			      tuplas[k]=indb[k];

			    this->Deformaciones->SetTuple(j,tuplas);
			    j++;
			 }		       
		    }
	       }
	  }
     }
 
    a=this->L-this->V;

   for (k=0;k<2;k++)
     this->EnE[k]->Initialize();
   this->Pg->Initialize();
/*   this->Puntos->Initialize();

   for (k=0;k<this->AlfaPot->GetNumberOfTuples();k++)
     this->Puntos->InsertNextPoint(this->AlfaPot->GetValue(k),this->L,this->V);
   for (k=0;k<this->Multiplicadores->GetNumberOfTuples();k++)
     this->Puntos->InsertNextPoint(this->Multiplicadores->GetComponent(k,0),this->Multiplicadores->GetComponent(k,1),this->Multiplicadores->GetComponent(k,2));*/
   
//   if (manera==0)
   this->O=this->O+1;
   this->P=manera;
   
   this->Rho->Initialize();
   this->Rho->SetScalarTypeToInt();
   this->Rho->SetNumberOfScalarComponents(this->NumEntidades);
   this->Rho->SetDimensions(this->DimMalla);
   
   RhoEnfr->Initialize();
   RhoEnfr->SetScalarTypeToInt();
   RhoEnfr->SetNumberOfScalarComponents(this->Ns);
   RhoEnfr->SetDimensions(this->DimMalla);

   //Nos interesa hacer zonas en cada componente de la malla.
   //Tambien nos interesa que los parametros pasen a ser un vtkStructuredPoints 
   //para hacer que reciban las magnitudes para cada zona.
   if (this->O<4)
     {
	this->Zonas->Initialize();
	this->Zonas->SetScalarTypeToInt();
	this->Zonas->SetNumberOfScalarComponents(3);
	this->Zonas->SetDimensions(this->DimZ[0]+1,this->DimZ[1]+1,this->DimZ[2]+1);
	this->Parametros->Initialize();
	this->Parametros->SetScalarTypeToFloat();
	this->Parametros->SetNumberOfScalarComponents(this->NumEntidades*this->L);
	this->Parametros->SetDimensions(this->DimZ);
     }
   else //Habria que diferenciar en el caso cardiaco en el caso de que decidamos coger las 17 regiones estandar.
     {
	this->Zonas->Initialize();
	this->Zonas->SetScalarTypeToInt();
	this->Zonas->SetNumberOfScalarComponents(3);
	this->Zonas->SetDimensions(this->DimZ[0]+1,this->DimZ[1]+1,this->DimZ[2]+1);
	this->Parametros->Initialize();
	this->Parametros->SetScalarTypeToFloat();
	this->Parametros->SetNumberOfScalarComponents(this->NumEntidades*this->L);
	this->Parametros->SetDimensions(this->DimZ);
     }
   
   for (k=0;k<2;k++)
     {
	this->EnE[k]->SetNumberOfComponents(1);
	this->EnE[k]->SetNumberOfTuples(this->K);
     }
 
   this->Pg->SetNumberOfComponents(1);
   this->Pg->SetNumberOfTuples(this->K);

   for (u=0;u<this->DimMalla[2];u++)
     {
	this->EnergiaNormalizada[u]->Initialize();
	this->Potencia[u]->Initialize();
	this->Potencia[u]->SetScalarTypeToFloat();
	this->EnergiaNormalizada[u]->SetScalarTypeToFloat();
	this->Potencia[u]->SetNumberOfScalarComponents(this->NumEntidades);
	this->EnergiaNormalizada[u]->SetNumberOfScalarComponents(this->NumEntidades*this->L);
	this->Potencia[u]->SetDimensions(this->K,this->DimMalla[0],this->DimMalla[1]);
	this->EnergiaNormalizada[u]->SetDimensions(this->K,this->DimMalla[0],this->DimMalla[1]);
     }

   if (this->P!=0 && !this->Leido)
     {
	//this->PuntParam->SetDataTypeToFloat();
	k=sprintf(nombre,"%s/vtk/Parametros.vtk",this->Ruta);
	Lector->SetFileName(nombre);
	Auxiliar=Lector->GetOutput();
	Auxiliar->Update();
	this->PuntParam=Auxiliar->GetPoints();
	kbis=this->PuntParam->GetNumberOfPoints()/this->NumEntidades;
        this->Leido=1;
	for (j=0;j<this->NumEntidades;j++)
	  {
	     for (n=0;n<this->V;n++)
	       paraux[n]=0.0;
	     for (k=0;k<kbis;k++)
	       {
		  this->PuntParam->GetPoint(k*this->NumEntidades+j,punto);
		  for (n=0;n<this->V;n++)
		    paraux[n]+=punto[n];
	       }
	     for (inda[2]=0;inda[2]<this->DimZ[2];inda[2]++)
	       {
		  for (inda[1]=0;inda[1]<this->DimZ[1];inda[1]++)
		    {
		       for (inda[0]=0;inda[0]<this->DimZ[0];inda[0]++)
			 {
			    for (n=0;n<this->V;n++)
			      this->Parametros->SetScalarComponentFromFloat(inda[0],inda[1],inda[2],j*this->L+a+n,sqrt(paraux[n]/(float)kbis));
			 }
		    }
	       }
	  }
	
	/*      for (k=0;k<this->NumEntidades*this->L;k++)
	paraux[k]=0.0;
      j=0;
       for (k=5;k<this->PuntParam->GetNumberOfPoints();k=k+6)
	   {
	      this->PuntParam->GetPoint(k,punto);
	      paraux[21]=paraux[21]+pow(punto[0],2);
	      this->PuntParam->GetPoint(k-1,punto);
	      paraux[20]=paraux[20]+pow(punto[2],2);
	      paraux[19]=paraux[19]+pow(punto[1],2);
	      paraux[14]=paraux[14]+pow(punto[0],2);
	      this->PuntParam->GetPoint(k-2,punto);
	      paraux[13]=paraux[13]+pow(punto[2],2);
	      paraux[12]=paraux[12]+pow(punto[1],2);
	      paraux[11]=paraux[11]+pow(punto[0],2);
	      this->PuntParam->GetPoint(k-3,punto);
	      paraux[10]=paraux[10]+pow(punto[0],2);
	      this->PuntParam->GetPoint(k-4,punto);
	      paraux[9]=paraux[9]+pow(punto[2],2);
	      paraux[8]=paraux[8]+pow(punto[1],2);
	      paraux[3]=paraux[3]+pow(punto[0],2);
	      this->PuntParam->GetPoint(k-5,punto);
	      paraux[2]=paraux[2]+pow(punto[2],2);
	      paraux[1]=paraux[1]+pow(punto[1],2);
	      paraux[0]=paraux[0]+pow(punto[0],2);
	      j=j+1;
	   }
      for (k=0;k<this->NumEntidades*this->L;k++)
	paraux[k]=paraux[k]/j;
      for (k=0;k<4;k++)
	{
      paraux[4+k]=paraux[k];
	   paraux[15+k]=paraux[11+k];
	}
      norma[0]=1+paraux[4]+paraux[5]+paraux[6]+paraux[7];
      norma[1]=1+paraux[15]+paraux[16]+paraux[17]+paraux[18];
      for (k=0;k<this->L;k++)
	paraux[k]=paraux[k]/norma[0];
      for (k=this->L;k<this->NumEntidades*this->L;k++)
	paraux[k]=paraux[k]/norma[1];

      for (k=0;k<this->NumEntidades*this->L;k++)
	this->Parametros->SetComponent(k,0,sqrt(paraux[k]));*/
	
     }
   
   for (inda[2]=0;inda[2]<this->DimZ[2]+1;inda[2]++)
     {
	for (inda[1]=0;inda[1]<this->DimZ[1]+1;inda[1]++)
	  {
	     for (inda[0]=0;inda[0]<this->DimZ[0]+1;inda[0]++)
	       {
		  for (j=0;j<3;j++)
		    this->Zonas->SetScalarComponentFromFloat(inda[0],inda[1],inda[2],j,(int)ceil(inda[j]*this->DimMalla[j]*1.0/this->DimZ[j]));
		  if (inda[0]!=this->DimZ[0] && inda[1]!=this->DimZ[1] && inda[2]!=this->DimZ[2])
		    {
//		       if (this->P==0)
//			 {
			    for (n=0;n<this->NumEntidades*this->L;n++)
			      this->Parametros->SetScalarComponentFromFloat(inda[0],inda[1],inda[2],n,this->Suma*sqrt(1.0/(float(this->L))));
/*			 }
		       else
			 {
			    for (k=0;k<this->NumEntidades;k++)
			      {
				 paraux[k]=0;
				 for (n=0;n<this->V;n++)
				   paraux[k]+=pow(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],k*this->L+a+n),2);
				 for (n=0;n<a;n++)
				   this->Parametros->SetScalarComponentFromFloat(inda[0],inda[1],inda[2],n+this->L*k,this->Suma*sqrt((1-paraux[k])/(float(a))));
			      }
			 }
		       for (z=0;z<this->NumEntidades*this->L-3;z+=3)
			 this->Puntos->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z+1),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z+2));
		       if ((this->NumEntidades*this->L)%3==0)
			 {
			    this->Puntos->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z+1),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z+2));
			    this->Puntos->InsertNextPoint(2,2,2);
			 }
		       else if ((this->NumEntidades*this->L)%3==1)
			 this->Puntos->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z),2,2);
		       else if ((this->NumEntidades*this->L)%3==2)
			 this->Puntos->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],z+1),2);*/
		    }
	       }
	  }
     }
   
   if (this->Super==1)
     {
	this->Beta=this->FactConv*log(4.0)/(log(this->Ns+5.0))/this->Suma;
	this->OptimizaParametros();
     }

   for (kbis=0;kbis<this->NumEntidades;kbis++)
     {
	for (inda[2]=0;inda[2]<this->DimZ[2];inda[2]++)
	  {
	     for (inda[1]=0;inda[1]<this->DimZ[1];inda[1]++)
	       {
		  for (inda[0]=0;inda[0]<this->DimZ[0];inda[0]++)
		    {
		       for (u=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],2);u<(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2]+1,2);u++)
			 {
			    for (t=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],1);t<(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1]+1,inda[2],1);t++)
			      {
				 for (j=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],0);j<(int)this->Zonas->GetScalarComponentAsFloat(inda[0]+1,inda[1],inda[2],0);j++)
				   {	
				      min=0;
				      for (n=0;n<this->V;n++)
					min+=this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],kbis*this->L+n+a)*this->LE[kbis][u]->GetScalarComponentAsFloat(0,j,t,n);
				      kmin=0;
				      for (k=1;k<this->K;k++)
					{
					   aux=0;
					   for (n=0;n<this->V;n++)
					     aux+=this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],kbis*this->L+n+a)*this->LE[kbis][u]->GetScalarComponentAsFloat(k,j,t,n);
					   if (aux<min)
					     {
						min=aux;
						kmin=k;
					     }
					}
				      this->Rho->SetScalarComponentFromFloat(j,t,u,kbis,kmin);
				      for (k=0;k<this->K;k++)
					{
					   if (k!=kmin)
					     this->CampoMedio[kbis]->SetScalarComponentFromFloat(j,t,u,k,0);
					   else
					     this->CampoMedio[kbis]->SetScalarComponentFromFloat(j,t,u,k,1);
					}				      
				   }
			      }			    
			 }
		    }
	       }
	  }
     }
   
   this->ConvierteRhoBordes();

   if (this->Super==0)
     {
	this->Beta=this->FactConv*log(4.0)/(log(this->Ns+5.0))/this->Suma;
	this->OptimizaParametros();
     }
   
//   this->Monitoriza();
   
   
   for (n=0;n<this->Ns;n+=this->Nsinter)
     {
	for (naux=0;naux<this->Nsinter;naux++)
	  {
	     this->Beta=this->FactConv*log(n+naux+5.0)/(log(this->Ns+5.0))/this->Suma;
	     for (u=0;u<this->DimMalla[0];u++)
	       {
		  for (j=0;j<this->DimMalla[1];j++)
		    RhoEnfr->SetScalarComponentFromFloat(u,j,0,n+naux,this->RhoBordes->GetScalarComponentAsFloat(u+this->Ord,j+this->Ord,this->Ord,0));
	       }
	     this->Ejecuta();
	     this->CampoMedio[0]->Update();
	     if (this->UsimagTool==0)
	       {		  
		  z=sprintf(nomb,"MeanField%d.vtk",n+naux);
		  EscribeRho->SetFileName(nomb);
		  EscribeRho->SetFileTypeToASCII();
		  EscribeRho->SetInput(this->CampoMedio[0]);
		  EscribeRho->Write();
	       }	     
	  }
	this->PuntObservacion->InsertNextPoint(n,n,n);
	this->PuntObservacion->InsertNextPoint(n,n,n);
	this->PuntObservacion->InsertNextPoint(n,n,n);
//	this->Monitoriza();
	if (this->Super==0)
	     this->OptimizaParametros();
     }
   //Caso renal:
   if (this->UsimagTool==0)
     {	
	z=sprintf(nomb,"%s.RhoEvol.vtk",this->Ruta);
	//Caso cardiaco:   
	//   z=sprintf(nomb,"%s/vtk/%dRhoEvol%d.vtk",this->Ruta,this->P,this->O);
	//   
	EscribeRho->SetFileName(nomb);
	EscribeRho->SetFileTypeToASCII();
	EscribeRho->SetInput(RhoEnfr);
	EscribeRho->Write();      
	
	z=sprintf(nomb,"%s/vtk/%dLENDO%d.vtk",this->Ruta,this->P,this->O);
	EscribeRho->SetFileName(nomb);
	EscribeRho->SetFileTypeToASCII();
	EscribeRho->SetInput(this->LE[0][0]);
	EscribeRho->Write();
   
	z=sprintf(nomb,"%s/vtk/%dLEPI%d.vtk",this->Ruta,this->P,this->O);
	EscribeRho->SetFileName(nomb);
	EscribeRho->SetFileTypeToASCII();
	EscribeRho->SetInput(this->LE[1][0]);
	EscribeRho->Write();
     }   
   
   this->ConvierteRho();
   if (this->P==0)
     {
/*	this->PuntParam->InsertNextPoint(this->Parametros->GetNumberOfTuples(),0,0);
	for (n=0;n<this->Parametros->GetNumberOfTuples()-3;n+=3)
	  this->PuntParam->InsertNextPoint(this->Parametros->GetComponent(n,0),this->Parametros->GetComponent(n+1,0),this->Parametros->GetComponent(n+2,0));
	if (this->Parametros->GetNumberOfTuples()%3==0)
	  {
	     this->PuntParam->InsertNextPoint(this->Parametros->GetComponent(n,0),this->Parametros->GetComponent(n+1,0),this->Parametros->GetComponent(n+2,0));
	     this->PuntParam->InsertNextPoint(2,2,2);
	  }
	else if (this->Parametros->GetNumberOfTuples()%3==1)
	  this->PuntParam->InsertNextPoint(this->Parametros->GetComponent(n,0),2,2);
	else if (this->Parametros->GetNumberOfTuples()%3==2)
	  this->PuntParam->InsertNextPoint(this->Parametros->GetComponent(n,0),this->Parametros->GetComponent(n+1,0),2);*/
	
	//En el caso renal habria que suspender la escritura hasta donde se indica. En el cardiaco no//
/*	this->PuntParam->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(0,0,0,a),this->Parametros->GetScalarComponentAsFloat(0,0,0,a+1),this->Parametros->GetScalarComponentAsFloat(0,0,0,a+2));
	this->PuntParam->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(0,0,0,this->L+a),this->Parametros->GetScalarComponentAsFloat(0,0,0,this->L+a+1),this->Parametros->GetScalarComponentAsFloat(0,0,0,this->L+a+2));
	Auxiliar->SetPoints(this->PuntParam);
	k=sprintf(nombre,"%s/vtk/Parametros.vtk",this->Ruta);
	Escritorc->SetFileName(nombre);
	Escritorc->SetFileTypeToASCII();
	Escritorc->SetInput(Auxiliar);
	Escritorc->Write();*/
	//Hasta aqui//
     }
   Poli->SetPoints(this->Puntos);
   if (this->UsimagTool==0)
     {	
	//Caso renal
	k=sprintf(nombre,"%s.Par.vtk",this->Ruta);
	//Caso cardiaco
//	k=sprintf(nombre,"%s/vtk/%dPar%d.vtk",this->Ruta,this->P,this->O);
	Escritora->SetFileName(nombre);
	Escritora->SetFileTypeToASCII();
	Escritora->SetInput(this->Parametros);
	Escritora->Write();
	
	//Caso renal
	z=sprintf(nombre,"%s.ParEvol.vtk",this->Ruta);
	//Caso cardiaco
//	z=sprintf(nombre,"%s/vtk/%dParEvol%d.vtk",this->Ruta,this->P,this->O);
	Escritorb->SetFileName(nombre);
	Escritorb->SetFileTypeToASCII();
	Escritorb->SetInput(Poli);
	Escritorb->Write();
     }     
   

//   this->PuntObservacion->InsertNextPoint(this->Puntos->GetNumberOfPoints(),10,11);
//   this->Monitoriza();
//for (k=0;k<this->Puntos->GetNumberOfPoints();k++)
//     {
//	this->PuntObservacion->InsertNextPoint(this->Puntos->GetPoint(k));
//	this->Monitoriza();
//     }
   //else
   //{
/*   Auxiliar->SetPoints(this->Puntos);
   k=sprintf(nombre,"MO%d.vtk",this->escribe);
   Escritor->SetFileName(nombre);
   Escritor->SetFileTypeToASCII();
   Escritor->SetInput(Auxiliar);
   Escritor->Write();
   //}
   this->escribe++;*/
}

void vtkOptimizaContorno::ConstruyeVecindario()
{
   bool existen;
   int di,el,au,num_punt,num_punu,i,j,k;
   int ele[4],inda[3],indb[3];
   int *punt,*punu;
   float aux;
   double punta[3],puntb[3],puntc[3];
   char nombre[200];
   
   vtkStructuredPoints *PPCategorizar=vtkStructuredPoints::New();
   vtkIntArray *LCategorizar=vtkIntArray::New();
   vtkIntArray *OCategorizar=vtkIntArray::New();
   vtkIntArray *PCategorizar=vtkIntArray::New();
   vtkCellArray *CliqAux=vtkCellArray::New();
   
   vtkPolyDataWriter *Escritor=vtkPolyDataWriter::New();
   vtkPolyData *Aux[3];
   for (i=0;i<this->DimMalla[2];i++)
     {
	this->EnergiaNormalizada[i]=vtkStructuredPoints::New();
	this->Potencia[i]=vtkStructuredPoints::New();
     }
   
   for (i=0;i<3;i++)
     Aux[i]=vtkPolyData::New();

   this->NumClNivel->Initialize();
   this->Vecindario->Initialize();
   this->Cliques->Initialize();
   this->Cliques->Reset();
   this->OrdenCliques->Initialize();
   this->OrdenCliques->Reset();
   this->CatCliques->Initialize();
   this->CatCliques->Reset();
   
   this->NumClNivel->SetNumberOfComponents(1);
   
   this->Ord=(int)floor(sqrt(this->OrdenMalla));
   
   this->Vecindario->InsertNextPoint(0,0,0);

   for (el=0;el<3;el++)
     inda[el]=0;
   for (di=0;di<this->DimensionalidadMalla;di++)
     inda[di]=1;
   for (indb[2]=-this->Ord*inda[2];indb[2]<=this->Ord*inda[2];indb[2]++)
     {
	for (indb[1]=-this->Ord*inda[1];indb[1]<=this->Ord*inda[1];indb[1]++)
	  {
	     for (indb[0]=-this->Ord*inda[0];indb[0]<=this->Ord*inda[0];indb[0]++)
	       {
		  if (this->Independencia==1)
		    {
		       if (((float)indb[0]!=0 && (float)indb[1]==0 && (float)indb[2]==0) || ((float)indb[0]==0 && (float)indb[1]!=0 && (float)indb[2]==0) || ((float)indb[0]==0 && (float)indb[1]==0 && (float)indb[2]!=0))
			 this->Vecindario->InsertNextPoint(indb[0],indb[1],indb[2]);
		    }
		  else
		    {  
		       aux=pow((float)indb[0],2)+pow((float)indb[1],2)+pow((float)indb[2],2);
		       if (aux!=0 && aux<=this->OrdenMalla)
			 au=this->Vecindario->InsertNextPoint(indb[0],indb[1],indb[2]);
		    }		       
	       }
	  }
     }
   
   for (i=0;i<3;i++)
     Aux[i]->SetPoints(this->Vecindario);
   
   this->Cliques->InitTraversal();
   this->Cliques->InsertNextCell(1);
   this->Cliques->InsertCellPoint(0);
   this->NumClNivel->InsertNextValue(1);
/*   if (this->Independencia==1)
     {
	for (i=1;i<=this->Ord;i++)
	  {
	     el=0;
	     for (k=0;k<this->DimensionalidadMalla;k++)
	       {
		  for (j=0;j<=i;j++)
		    {
		       this->Cliques->InsertNextCell(i+1);
		       this->Cliques->InsertCellPoint(0);
		       for (di=0;di<i;di++)
			 this->Cliques->InsertCellPoint(this->Ord+1-i+j+di+2*this->Ord*k);
		    }
		  el+=i+1;
	       }
	     this->NumClNivel->InsertNextValue(el);
	  }
     }*/
//   else
//     {
	existen=0;
	for (au=1;au<this->Vecindario->GetNumberOfPoints();au++)
	  {
	     this->Cliques->InsertNextCell(2);
	     this->Cliques->InsertCellPoint(0);
	     this->Cliques->InsertCellPoint(au);
	     existen=1;
	  }
	this->NumClNivel->InsertNextValue(this->Vecindario->GetNumberOfPoints()-1);
	au=1;
	while (existen && au<2)
	  {
		 CliqAux->Initialize();
		 CliqAux->Reset();
		 CliqAux->DeepCopy(this->Cliques);
	     existen=0;
	     ele[0]=0;
	     ele[2]=0;
	     for (el=0;el<au;el++)
	       ele[0]+=this->NumClNivel->GetValue(el);
	     ele[1]=ele[0]+this->NumClNivel->GetValue(au);
	     au++;
	     LCategorizar->Initialize();
	     LCategorizar->SetNumberOfComponents(1);
	     LCategorizar->SetNumberOfTuples(ele[1]-ele[0]);
	     for (el=ele[0];el<ele[1];el++)
	       {
		  this->Cliques->GetNextCell(num_punt,punt);
		  LCategorizar->SetValue(el-ele[0],this->Cliques->GetTraversalLocation());
	       }
	     for (el=ele[0];el<ele[1];el++)
	       {
		  CliqAux->GetCell(LCategorizar->GetValue(el-ele[0]),num_punt,punt);
		  this->Vecindario->GetPoint(punt[num_punt-1],punta);
		  for (di=punt[num_punt-1]+1;di<this->Vecindario->GetNumberOfPoints();di++) //Ordenacion estrictamente creciente.
		    {
		       this->Vecindario->GetPoint(di,puntb);
		       aux=0;
		       for(i=0;i<3;i++)
			 {
			    puntb[i]=pow(puntb[i]-punta[i],2);
			    aux+=puntb[i];
			 }
		       if (aux<=this->OrdenMalla) //Condicion de que sean vecinos.
			 {
			    for (j=ele[0];j<ele[1];j++)
			      {
				 CliqAux->GetCell(LCategorizar->GetValue(j-ele[0]),num_punu,punu);
				 //Esto parece que no esta bien; lo de arriba y lo de abajo.
				 if (punu[num_punu-2]==punt[num_punt-2] && punu[num_punu-1]==di) //Condicion de fuentes comunes con el nivel anterior.
				   {
				      if (this->Independencia==0 || (this->Independencia==1 && ((puntb[0]!=0 && puntb[1]==0 && puntb[2]==0) || (puntb[0]==0 && puntb[1]!=0 && puntb[2]==0) || (puntb[0]==0 && puntb[1]==0 && puntb[2]!=0))))
					{					   
					   existen=1;
					   this->Cliques->InsertNextCell(num_punt+1);
					   for (i=0;i<num_punt;i++)
					     this->Cliques->InsertCellPoint(punt[i]);
					   this->Cliques->InsertCellPoint(di);
					   ele[2]++;
					   break;
					}				      
				   }
			      }
			 }
		    }
	       }
	     this->NumClNivel->InsertNextValue(ele[2]);
	  }
//     }
   
   for (i=0;i<this->NumClNivel->GetNumberOfTuples();i++)
     this->NumClNivel->SetValue(i,this->NumClNivel->GetValue(i)/(i+1));  //De esta manera tengo el numero de categorias de cliques por nivel.
   
   this->Multiplicadores->Initialize();
   this->Multiplicadores->SetNumberOfComponents(1);

   //Se ordenan los cliques convenientemente.
   this->Cliques->InitTraversal();
   for (i=0;i<this->Cliques->GetNumberOfCells();i++)
     {
	this->Cliques->GetNextCell(num_punt,punt);
	this->Multiplicadores->InsertNextValue(0);
	LCategorizar->Initialize();
	LCategorizar->SetNumberOfComponents(1);
	LCategorizar->SetNumberOfTuples(num_punt);
	for (j=0;j<num_punt;j++)
	  LCategorizar->SetValue(j,j);
	
	for (j=0;j<num_punt;j++)
	  {
	     this->Vecindario->GetPoint(punt[LCategorizar->GetValue(j)],punta);
	     for (k=j+1;k<num_punt;k++)
	       {
		  this->Vecindario->GetPoint(punt[LCategorizar->GetValue(k)],puntb);
		  if (puntb[0]<punta[0] || (puntb[0]==punta[0] && (puntb[1]<punta[1] || (puntb[1]==punta[1] && puntb[2]<punta[2]))))
		    {
		       el=LCategorizar->GetValue(k);
		       LCategorizar->SetValue(k,LCategorizar->GetValue(j));
		       LCategorizar->SetValue(j,el);
		       for (di=0;di<3;di++)
			 punta[di]=puntb[di];
		    }
	       }
	  }
	if (num_punt==3)
	  {
	     di=0;
	     for (j=0;j<num_punt;j++)
	       {
		  this->Vecindario->GetPoint(punt[LCategorizar->GetValue(j)],punta);
		  for (k=j+1;k<num_punt;k++)
		    {
		       this->Vecindario->GetPoint(punt[LCategorizar->GetValue(k)],puntb);
		       puntc[di]=sqrt(pow(puntb[0]-punta[0],2)+pow(puntb[1]-punta[1],2)+pow(puntb[2]-punta[2],2));
		       di++;
		    }
	       }
	     OCategorizar->Initialize();
	     OCategorizar->SetNumberOfComponents(1);
	     OCategorizar->SetNumberOfTuples(di);
	     for (j=0;j<di;j++)
	       OCategorizar->SetValue(j,j);
	     for (j=0;j<di-1;j++)
	       {
		  for (k=j+1;k<di;k++)
		    {
		       if (puntc[k]<puntc[j])
			 {
			    aux=puntc[k];
			    puntc[k]=puntc[j];
			    puntc[j]=aux;
			    el=OCategorizar->GetValue(k);
			    OCategorizar->SetValue(k,OCategorizar->GetValue(j));
			    OCategorizar->SetValue(j,el);
			 }
		    }
	       }
	     el=OCategorizar->GetValue(0);
	     for (k=num_punt-1;k>0;k--)
	       {
		  el=el-k;
		  if (el<0)
		    {
		       ele[0]=num_punt-1-k;
		       ele[1]=num_punt+el;
		       break;
		    }
	       }
	     el=OCategorizar->GetValue(1);
	     for (k=num_punt-1;k>0;k--)
	       {
		  el=el-k;
		  if (el<0)
		    {
		       ele[2]=num_punt-1-k;
		       ele[3]=num_punt+el;
		       break;
		    }
	       }
	     if (ele[0]==ele[2] || ele[0]==ele[3])
	       {
		  if (ele[0]==ele[2])
		    {
		       OCategorizar->SetValue(0,LCategorizar->GetValue(ele[1]));
		       OCategorizar->SetValue(1,LCategorizar->GetValue(ele[0]));
		       OCategorizar->SetValue(2,LCategorizar->GetValue(ele[3]));
		    }
		  else
		    {
		       OCategorizar->SetValue(0,LCategorizar->GetValue(ele[1]));
		       OCategorizar->SetValue(1,LCategorizar->GetValue(ele[0]));
		       OCategorizar->SetValue(2,LCategorizar->GetValue(ele[2]));
		    }
	       }
	     else
	       {
		  if (ele[1]==ele[2])
		    {
		       OCategorizar->SetValue(0,LCategorizar->GetValue(ele[0]));
		       OCategorizar->SetValue(1,LCategorizar->GetValue(ele[1]));
		       OCategorizar->SetValue(2,LCategorizar->GetValue(ele[3]));
		    }
		  else
		    {
		       OCategorizar->SetValue(0,LCategorizar->GetValue(ele[0]));
		       OCategorizar->SetValue(1,LCategorizar->GetValue(ele[1]));
		       OCategorizar->SetValue(2,LCategorizar->GetValue(ele[2]));
		    }
	       }
	     for (k=0;k<num_punt;k++)
	       LCategorizar->SetValue(k,OCategorizar->GetValue(k));
	  }
	this->OrdenCliques->InsertNextCell(num_punt);
	for (j=0;j<num_punt;j++)
	  this->OrdenCliques->InsertCellPoint(LCategorizar->GetValue(j));
	if (num_punt==1)
	  this->Multiplicadores->InsertNextValue(1);
	else if (num_punt==2)
	  {
	     this->Vecindario->GetPoint(punt[LCategorizar->GetValue(0)],punta);
	     this->Vecindario->GetPoint(punt[LCategorizar->GetValue(1)],puntb);
	     aux=sqrt(pow(puntb[0]-punta[0],2)+pow(puntb[1]-punta[1],2)+pow(puntb[2]-punta[2],2));
	     this->Multiplicadores->InsertNextValue(-1.0/aux);
	     this->Multiplicadores->InsertNextValue(1.0/aux);
	  }
	else
	  {
	     this->Vecindario->GetPoint(punt[LCategorizar->GetValue(0)],punta);
	     this->Vecindario->GetPoint(punt[LCategorizar->GetValue(1)],puntb);
	     puntc[0]=sqrt(pow(puntb[0]-punta[0],2)+pow(puntb[1]-punta[1],2)+pow(puntb[2]-punta[2],2));	     
	     this->Multiplicadores->InsertNextValue(1.0/puntc[0]);
	     this->Vecindario->GetPoint(punt[LCategorizar->GetValue(1)],punta);
	     this->Vecindario->GetPoint(punt[LCategorizar->GetValue(2)],puntb);
	     puntc[1]=sqrt(pow(puntb[0]-punta[0],2)+pow(puntb[1]-punta[1],2)+pow(puntb[2]-punta[2],2));
	     this->Multiplicadores->InsertNextValue(-1.0/puntc[0]-1.0/puntc[1]);
	     this->Multiplicadores->InsertNextValue(1.0/puntc[1]);
	  }
     }

/*   if (this->Independencia==1)
     {	
	//Clasificacion de los cliques en categorias.
	this->Cliques->InitTraversal();
	for (i=0;i<this->NumClNivel->GetNumberOfTuples();i++)
	  {
	     for (j=0;j<this->NumClNivel->GetValue(i);j++)
	       {
		  this->CatCliques->InsertNextCell(i+1);
		  for (k=0;k<i+1;k++)
		    {
		       this->CatCliques->InsertCellPoint(this->Cliques->GetTraversalLocation());
		       this->Cliques->GetNextCell(num_punt,punt);
		    }
	       }
	  }
     }*/
//   else
//     {
	ele[0]=0;
	this->Cliques->InitTraversal();
	this->OrdenCliques->InitTraversal();
	//Clasificacion de los cliques en categorias.
	for (i=0;i<this->NumClNivel->GetNumberOfTuples();i++)
	  {
	     ele[1]=ele[0]+this->NumClNivel->GetValue(i)*(i+1);
	     PPCategorizar->Initialize();
	     LCategorizar->Initialize();
	     PCategorizar->Initialize();
	     PPCategorizar->SetScalarTypeToInt();
	     LCategorizar->SetNumberOfComponents(1);
	     LCategorizar->SetNumberOfTuples(ele[1]-ele[0]);
	     PCategorizar->SetNumberOfComponents(1);
	     PCategorizar->SetNumberOfTuples(ele[1]-ele[0]);
	     PPCategorizar->SetNumberOfScalarComponents(1);
	     PPCategorizar->SetDimensions(ele[1]-ele[0],i,3);
	     for (di=ele[0];di<ele[1];di++)
	       {
		  LCategorizar->SetValue(di-ele[0],0);
		  PCategorizar->SetValue(di-ele[0],this->Cliques->GetTraversalLocation());
		  this->Cliques->GetNextCell(num_punt,punt);
		  this->OrdenCliques->GetNextCell(num_punu,punu);
		  for (j=0;j<num_punt-1;j++)
		    {
		       this->Vecindario->GetPoint(punt[punu[j]],punta);
		       this->Vecindario->GetPoint(punt[punu[j+1]],puntb);
		       for (ele[2]=0;ele[2]<3;ele[2]++)
			 {
			    puntb[ele[2]]=puntb[ele[2]]-punta[ele[2]];
			    PPCategorizar->SetScalarComponentFromFloat(di-ele[0],j,ele[2],0,puntb[ele[2]]);
			 }
		    }
		  //Para hacerlo hay que restar todos los puntos de la celda original entre si y ver cuales cuadran emparejados.
	       }
	     for (di=ele[0];di<ele[1];di++)
	       {
		  if (LCategorizar->GetValue(di-ele[0])==0)
		    {
		       this->CatCliques->InsertNextCell(i+1);
		       this->CatCliques->InsertCellPoint(PCategorizar->GetValue(di-ele[0]));
		       LCategorizar->SetValue(di-ele[0],1);
		       OCategorizar->Initialize();
		       OCategorizar->SetNumberOfComponents(1);
		       OCategorizar->SetNumberOfTuples(num_punt-1);
		       for (el=di+1;el<ele[1];el++)
			 {
			    k=0;
			    for (j=0;j<num_punt-1;j++)
			      OCategorizar->SetValue(j,0);
			    for (j=0;j<num_punt-1;j++)
			      {
				 for (ele[2]=0;ele[2]<3;ele[2]++)
				   {
				      punta[ele[2]]=PPCategorizar->GetScalarComponentAsFloat(di-ele[0],j,ele[2],0);
				      puntb[ele[2]]=PPCategorizar->GetScalarComponentAsFloat(el-ele[0],j,ele[2],0);
				   }
				 if (punta[0]==puntb[0] && punta[1]==puntb[1] && punta[2]==puntb[2] && OCategorizar->GetValue(k)==0)
				   {
				      OCategorizar->SetValue(k,1);
				      k++;
				   }	
			      }
			    if (k==num_punt-1)
			      {
				 this->CatCliques->InsertCellPoint(PCategorizar->GetValue(el-ele[0]));
				 LCategorizar->SetValue(el-ele[0],1);
			      }
			 }
		    }
	       }
	     ele[0]=ele[1];
	  }
//     }
   
   this->L=this->V;
   if (this->NumEntidades==1)
     {
	for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
	  this->L+=this->NumClNivel->GetValue(i);
     }
   else
     {
	for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
	  this->L+=2*this->NumClNivel->GetValue(i);
     }
   
   this->AlfaPot->SetNumberOfComponents(1);
   this->BetaPot->SetNumberOfComponents(1);
   this->AlfaPot->SetNumberOfTuples(this->L);
   di=this->L-this->V;
   this->BetaPot->SetNumberOfTuples(di);
   
   for (i=0;i<di;i++)
     this->BetaPot->SetValue(i,.85);
   //Para el caso renal, este valor funciona bien en .9 para un salto=3.
   
   /////TENGO UN PROBLEMA CON TODO ESTO, PUES CAMBIE LO QUE CAMBIE, NO PARECE AFECTAR AL RESULTADO FINAL; HAY QUE BUSCAR VARIACIONES DRASTICAS EN EL CASO TEMPORAL.
   if (this->NumEntidades==1)
     {
	el=0;
	for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
	  {
	     for (j=0;j<this->NumClNivel->GetValue(i);j++)
	       {
		  this->AlfaPot->SetValue(el,1.0/(pow(2.0,i)));
		  el++;
	       }
	  }
     }
   else
     {
	el=0;
	for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
	  {
	     for (j=0;j<this->NumClNivel->GetValue(i);j++)
	       {
		  this->AlfaPot->SetValue(el,1.0/(pow(2.0,i)));
		  this->AlfaPot->SetValue(di/2+el,1.0/(2.0*pow(2.0,i)));
		  el++;
	       }
	  }
     }
   for (i=di;i<this->AlfaPot->GetNumberOfTuples();i++)
     this->AlfaPot->SetValue(i,1.0);

/*   this->ParamVeros->SetNumberOfComponents(2);
   this->ParamVeres->SetNumberOfTuples(this->V);
   for (i=0;i<this->NumEntidades;i++)
     {
	for (j=0;j<this->V;j++)
	  this->ParamVeros->SetComponent(j,i,0);
     }*/
   
   
   //Se construyen los multiplicadores del metodo de diferencias finitas
 /*  this->Multiplicadores->SetNumberOfComponents(this->NumClNivel->GetNumberOfTuples());
   this->Multiplicadores->SetNumberOfTuples(this->NumClNivel->GetNumberOfTuples());
   
   //Habria que revisar este metodo de construccion de los multiplicadores, sobre todo para casos en que el
   //NumClNivel->Tuples sea mayor que 3.
   for (i=0;i<this->NumClNivel->GetNumberOfTuples();i++)
     {
	for (j=0;j<this->NumClNivel->GetNumberOfTuples();j++)
	  {
	     this->Multiplicadores->SetComponent(i,j,0);
	  }
     }
   this->Multiplicadores->SetComponent(0,0,1);
   for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
     {
	this->Multiplicadores->SetComponent(i,0,-this->Multiplicadores->GetComponent(i-1,0));
	for (j=1;j<=i;j++)
	  this->Multiplicadores->SetComponent(i,j,this->Multiplicadores->GetComponent(i-1,j-1)-this->Multiplicadores->GetComponent(i-1,j));
     }*/
   
//   if (this->escribe==1)
 //    {

   if (this->UsimagTool==0)
     {	
	//Caso renal:
	i=sprintf(nombre,"%s.CliqEl.vtk",this->Ruta);
	Aux[0]->SetVerts(this->Cliques);
	Escritor->SetFileName(nombre);
	Escritor->SetFileTypeToASCII();
	Escritor->SetInput(Aux[0]);
	Escritor->Write();
   
   
   //Aux[1]->SetVerts(this->OrdenCliques);
   //Escritor->SetFileName("ConstruyeCliquesb.vtk");
   //Escritor->SetFileTypeToASCII();
   //Escritor->SetInput(Aux[1]);
   //Escritor->Write();
   
	i=sprintf(nombre,"%s.CliqCat.vtk",this->Ruta);
	Aux[2]->SetVerts(this->CatCliques);
   	Escritor->SetFileName(nombre);
	Escritor->SetFileTypeToASCII();
	Escritor->SetInput(Aux[2]);
	Escritor->Write();
     }
   
//	this->escribe=0;
//     }
   
}

void vtkOptimizaContorno::ConstruyeMetricaCliques()
{
   int i,num_p,*num_pu,ORD=0;

   this->OrdenCliques->InitTraversal();
   this->OrdenCliques->GetNextCell(num_p,num_pu);
   for (i=0;i<this->OrdenCliques->GetNumberOfCells()-1;i++)
     {
	this->OrdenCliques->GetNextCell(num_p,num_pu);
	ORD+=num_p;
     }

   this->MetricaCliques->Initialize();
   this->MetricaCliques->SetScalarTypeToFloat();
   this->MetricaCliques->SetNumberOfScalarComponents(ORD);
   this->MetricaCliques->SetDimensions(this->DimMalla);
   
   //AQUI.
   //El problema es que necesitamos un campo de posiciones de la malla para todo esto.
   //La manera de definirlo es generarlo en vtkFuncionVerosimilitud O MAS BIEN EN VTKPLANTILLAAJUSTADA 
   //(DONDE SE GENERA EL RHO NULO QUE NOS INTERESA).
   //Puede ser un a√±adido a la estructura LR.
   //Bastaria con el campo para k=K/2 o necesitamos todos los k?.
   //Parece que bastaria para k=K/2. El problema seria la adimensionalidad en la coordenada radial, para lo que 
   //se puede usar la estimacion de parametros por zonas.
   
}

void vtkOptimizaContorno::OptimizaParametros()
{  
   bool entrada;
   int v,j,s,k,KE,jbis,t,u,inda[3],indb[3];
   float r;
   double temperatura;
   double candidato,rango[2],probact,probsig,numeract,numersig,denomact,denomsig,suma;
 
   
   vtkPolyDataWriter *Escritor=vtkPolyDataWriter::New();
   vtkPolyData *Poli=vtkPolyData::New();
   vtkPoints *Puntoss=vtkPoints::New();
   
   this->CalculaEnergiaNormalizada();
   KE=0;
   
   for (v=0;v<this->Nsparam;v++)
     {
	temperatura=16*log(v+2.0)/this->Suma;//(log(this->Nsparam+5.0));
	for (inda[2]=0;inda[2]<this->DimZ[2];inda[2]++)
	  {
	     for (inda[1]=0;inda[1]<this->DimZ[1];inda[1]++)
	       {
		  for (inda[0]=0;inda[0]<this->DimZ[0];inda[0]++)
		    {
		       for (j=0;j<this->NumEntidades*this->L;j++)
			 {
			    //2 corresponde a lo que en el algoritmo se llama 
			    //gamma; por el proyecto de RocÌo se entiende que
			    //Marcos ha desarrollado un mÈtodo que optimiza 
			    //esta elecciÛn est·tica del par·metro haciÈndolo 
			    //variar durante la optimizaciÛn seg˙n unas determinadas
			    //recetas; por otro lado, en el libro de Redes Neuronales
			    //tambiÈn se mostraban maneras de trabajar con 
			    //intervalos de descenso del gradiente adaptativos (creo).
			    //MIRAR TODO ESTO CUANDO TENGAMOS EL MARCO 
			    //BASE CONSTRUIDO
			    /*if (this->Super==1)
			     {
			     rango[0]=this->Parametros->GetValue(j)-(10*this->Parametros->GetValue(j)+1)/100.0;
			     rango[1]=this->Parametros->GetValue(j)+(10*this->Parametros->GetValue(j)+1)/100.0;
			     }
			     else
			     {*/
			    rango[0]=this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j)-(10*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j)+1)/10.0;
			    rango[1]=this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j)+(10*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j)+1)/10.0;
			    //}
			    r=(float)(1.0*rand()/(RAND_MAX+1.0));
			    candidato=r*(rango[1]-rango[0])+rango[0];
			    if (candidato<0)
			      candidato=0;
			    if (candidato>this->Suma)
			      candidato=this->Suma;
			    probact=0;
			    probsig=0;
			    numeract=0;
			    numersig=0;
			    suma=-this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j);
			    if (j<this->L)
			      {
				 for (jbis=0;jbis<this->L;jbis++)
				   suma=suma+this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],jbis)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],jbis);
				 suma=suma+candidato*candidato;
				 suma=sqrt(suma);
			      }
			    else
			      {
				 for (jbis=this->L;jbis<2*this->L;jbis++)
				   suma=suma+this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],jbis)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],jbis);
				 suma=suma+candidato*candidato;
				 suma=sqrt(suma);
			      }
			    for (u=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],2);u<(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2]+1,2);u++)
			      {
				 for (t=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],1);t<(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1]+1,inda[2],1);t++)
				   {
				      for (s=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],0);s<(int)this->Zonas->GetScalarComponentAsFloat(inda[0]+1,inda[1],inda[2],0);s++)
					{
					   denomact=0;
					   denomsig=0;
					   if (j>=this->L)
					     {
						if (this->Super==0)
						  KE=(int)this->RhoBordes->GetScalarComponentAsFloat(s+this->Ord,t+this->Ord,u+this->Ord,1);
						else
						  KE=(int)this->RhoMarMed->GetScalarComponentAsFloat(s+this->Ord,t,u,1); 
						
						numeract=numeract+this->Potencia[u]->GetScalarComponentAsFloat(KE,s,t,1);
						numersig=numersig+this->Potencia[u]->GetScalarComponentAsFloat(KE,s,t,1)+this->EnergiaNormalizada[u]->GetScalarComponentAsFloat(KE,s,t,j)*(candidato-this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j));
						for (k=0;k<this->K;k++)
						  {     
						     denomact=denomact+exp(this->Potencia[u]->GetScalarComponentAsFloat(k,s,t,1));
						     denomsig=denomsig+exp((this->Potencia[u]->GetScalarComponentAsFloat(k,s,t,1)+this->EnergiaNormalizada[u]->GetScalarComponentAsFloat(k,s,t,j)*(candidato-this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j)))/suma);
						  }
					     }
					   else
					     {
						if (this->Super==0)
						  KE=(int)this->RhoBordes->GetScalarComponentAsFloat(s+this->Ord,t+this->Ord,u+this->Ord,0);
						else
						  KE=(int)this->RhoMarMed->GetScalarComponentAsFloat(s+this->Ord,t,u,0);
						
						numeract=numeract+this->Potencia[u]->GetScalarComponentAsFloat(KE,s,t,0);
						numersig=numersig+this->Potencia[u]->GetScalarComponentAsFloat(KE,s,t,0)+this->EnergiaNormalizada[u]->GetScalarComponentAsFloat(KE,s,t,j)*(candidato-this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j));
						for (k=0;k<this->K;k++)
						  { 
						     denomact=denomact+exp(this->Potencia[u]->GetScalarComponentAsFloat(k,s,t,0));
						     denomsig=denomsig+exp((this->Potencia[u]->GetScalarComponentAsFloat(k,s,t,0)+this->EnergiaNormalizada[u]->GetScalarComponentAsFloat(k,s,t,j)*(candidato-this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j)))/suma);
						  }
					     }	
					   probact=probact+log(denomact);
					   probsig=probsig+log(denomsig);
					}
				   }
			      }			    
			    numersig=numersig/suma;
			    this->PuntObservacion->InsertNextPoint(numeract,probact,numersig);
			    this->PuntObservacion->InsertNextPoint(probsig,denomact,denomsig);
			    probact=exp(-(numeract-probact-numersig+probsig)*temperatura*this->Beta);
			    this->PuntObservacion->InsertNextPoint(probact,v,j);
			    entrada=1;
			    if (probact<1)
			      {
				 r=(float)(1.0*rand()/(RAND_MAX+1.0));
				 if (r>probact)
				   entrada=0;
			      }
			    if (entrada==1)
			      {
				 this->Parametros->SetScalarComponentFromFloat(inda[0],inda[1],inda[2],j,candidato);
				 if (j<this->L)
				   {
				      for (jbis=0;jbis<this->L;jbis++)
					this->Parametros->SetScalarComponentFromFloat(inda[0],inda[1],inda[2],jbis,this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],jbis)/suma);
				   }
				 else
				   {
				      for (jbis=this->L;jbis<2*this->L;jbis++)
					this->Parametros->SetScalarComponentFromFloat(inda[0],inda[1],inda[2],jbis,this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],jbis)/suma);				 
				   }
				 for (indb[2]=0;indb[2]<this->DimZ[2];indb[2]++)
				   {
				      for (indb[1]=0;indb[1]<this->DimZ[1];indb[1]++)
					{
					   for (indb[0]=0;indb[0]<this->DimZ[0];indb[0]++)
					     {
						for (u=(int)this->Zonas->GetScalarComponentAsFloat(indb[0],indb[1],indb[2],2);u<(int)this->Zonas->GetScalarComponentAsFloat(indb[0],indb[1],indb[2]+1,2);u++)
						  {
						     for (t=(int)this->Zonas->GetScalarComponentAsFloat(indb[0],indb[1],indb[2],1);t<(int)this->Zonas->GetScalarComponentAsFloat(indb[0],indb[1]+1,indb[2],1);t++)
						       {
							  for (s=(int)this->Zonas->GetScalarComponentAsFloat(indb[0],indb[1],indb[2],0);s<(int)this->Zonas->GetScalarComponentAsFloat(indb[0]+1,indb[1],indb[2],0);s++)
							    {
							       for (k=0;k<this->K;k++)
								 {
								    r=0;
								    if (j<this->L)
								      {
									 for (jbis=0;jbis<this->L;jbis++)
									   r=r+this->EnergiaNormalizada[u]->GetScalarComponentAsFloat(k,s,t,jbis)*this->Parametros->GetScalarComponentAsFloat(indb[0],indb[1],indb[2],jbis);
									 this->Potencia[u]->SetScalarComponentFromFloat(k,s,t,0,r);
								      }
								    else
								      {
									 for (jbis=this->L;jbis<2*this->L;jbis++)
									   r=r+this->EnergiaNormalizada[u]->GetScalarComponentAsFloat(k,s,t,jbis)*this->Parametros->GetScalarComponentAsFloat(indb[0],indb[1],indb[2],jbis);
									 this->Potencia[u]->SetScalarComponentFromFloat(k,s,t,1,r);
								      }
								 }
							    }
						       }  
						  }
					     }
					}
				   }
			      }
			 }
		       for (j=0;j<this->L*this->NumEntidades-3;j+=3)
			 this->Puntos->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j+1),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j+2));
		       if ((this->L*this->NumEntidades)%3==0)
			 {
			    this->Puntos->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j+1),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j+2));
			    this->Puntos->InsertNextPoint(2,2,2);
			 }
		       else if ((this->L*this->NumEntidades)%3==1)
			 this->Puntos->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j),2,2);
		       else if ((this->L*this->NumEntidades)%3==2)
			 this->Puntos->InsertNextPoint(this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j),this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j+1),2);
		    }
	       }
	  }
     }
   
   //   this->escribe++;
   Puntoss->Delete();
   Poli->Delete();
   Escritor->Delete();
}

void vtkOptimizaContorno::CalculaEnergiaNormalizada()
{
   //s==N˙mero de rayo
   //k==N˙mero de etiqueta
   int s,k,j,l,m,n,o,a,t,i,p,u,num_p,num_q,num_r, *num_pu, *num_qu, *num_ru,inda[3];
   float pot, aux, maxE;
   float *pot_aux;
   double punt[3];
   vtkFloatArray *Vecinos=vtkFloatArray::New();
   pot_aux = new float[this->DimensionalidadEstados];

   Vecinos->SetNumberOfComponents(this->NumEntidades*this->DimensionalidadEstados);
   Vecinos->SetNumberOfTuples(this->Vecindario->GetNumberOfPoints());
   a=this->L-this->V;
   
   for (inda[2]=0;inda[2]<this->DimZ[2];inda[2]++)
     {
	for (inda[1]=0;inda[1]<this->DimZ[1];inda[1]++)
	  {
	     for (inda[0]=0;inda[0]<this->DimZ[0];inda[0]++)
	       {
		  for (u=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],2);u<(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2]+1,2);u++)
		    {
		       for (t=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],1);t<(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1]+1,inda[2],1);t++)
			 {
			    for (s=(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],0);s<(int)this->Zonas->GetScalarComponentAsFloat(inda[0]+1,inda[1],inda[2],0);s++)
			      {
				   {
				      if (this->Super==0)
					{
					   for (j=0;j<this->NumEntidades;j++)
					     {
						for (l=0;l<this->DimensionalidadEstados;l++)
						  {
						     for (i=1;i<this->Vecindario->GetNumberOfPoints();i++)
						       {
							  this->Vecindario->GetPoint(i,punt);
							  pot=this->Deformaciones->GetComponent((int)this->RhoBordes->GetScalarComponentAsFloat(s+this->Ord+(int)punt[0],t+this->Ord+(int)punt[1],u+this->Ord+(int)punt[2],j),l);
							  Vecinos->SetComponent(i,j+l*this->DimensionalidadEstados,pot);
						       }
						  }
					     }
					}
				      else
					{
					   for (j=0;j<this->NumEntidades;j++)
					     {
						for (l=0;l<this->DimensionalidadEstados;l++)
						  {
						     for (i=1;i<this->Vecindario->GetNumberOfPoints();i++)
						       {
							  this->Vecindario->GetPoint(i,punt);
							  pot=this->Deformaciones->GetComponent((int)this->RhoMarMed->GetScalarComponentAsFloat(s+this->Ord+(int)punt[0],t+(int)punt[1],u+(int)punt[2],j),l);
							  Vecinos->SetComponent(i,j+l*this->DimensionalidadEstados,pot);
						       }  
						  }
					     }
					}
				      for (j=0;j<this->NumEntidades;j++)
					{
					   for (l=0;l<this->DimensionalidadEstados;l++)
					     {
						for (i=0;i<this->NumEntidades;i++)
						  {
						     if (i!=j)
						       {
							  if (this->Super==0)
							    pot=this->Deformaciones->GetComponent((int)this->RhoBordes->GetScalarComponentAsFloat(s+this->Ord,t+this->Ord,u+this->Ord,i),l);
							  else
							    pot=this->Deformaciones->GetComponent((int)this->RhoMarMed->GetScalarComponentAsFloat(s+this->Ord,t,u,i),l);
							  Vecinos->SetComponent(0,i+l*this->DimensionalidadEstados,pot);
						       }
						  }
					     }
					   for (k=0;k<this->K;k++)
					     {
						for (l=0;l<this->DimensionalidadEstados;l++)
						  {
						     pot=this->Deformaciones->GetComponent(k,l);
						     Vecinos->SetComponent(0,j+l*this->DimensionalidadEstados,pot);
						  }
						//Hay que hacer bien el potencial una vez que estan bien los vecinos.
						p=0;
						this->CatCliques->InitTraversal();
						this->CatCliques->GetNextCell(num_p,num_pu);
						for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
						  {
						     for (l=0;l<this->NumClNivel->GetValue(i);l++)
						       {
							  p++;
							  this->CatCliques->GetNextCell(num_p,num_pu);
							  pot=0;
							  for (m=0;m<num_p;m++)
							    {
							       this->Cliques->GetCell(num_pu[m],num_q,num_qu);
							       this->OrdenCliques->GetCell(num_pu[m],num_r,num_ru);
							       for (o=0;o<this->DimensionalidadEstados;o++)
								 {
								    pot_aux[o]=0;
								    for (n=0;n<num_r;n++)
								      pot_aux[o]+=Vecinos->GetComponent(num_qu[num_ru[n]],j+o*this->DimensionalidadEstados)*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
								 }
							       pot+=this->Potencial(pot_aux,p-1);
							    }
							  this->EnergiaNormalizada[u]->SetScalarComponentFromFloat(k,s,t,j*this->L+p-1,pot);
						       }
						  }
						if (this->NumEntidades>1)
						  {
						     p=0;
						     this->CatCliques->InitTraversal();
						     this->CatCliques->GetNextCell(num_p,num_pu);
						     for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
						       {
							  for (l=0;l<this->NumClNivel->GetValue(i);l++)
							    {
							       p++;
							       this->CatCliques->GetNextCell(num_p,num_pu);
							       pot=0;
							       for (m=0;m<num_p;m++)
								 {
								    this->Cliques->GetCell(num_pu[m],num_q,num_qu);
								    this->OrdenCliques->GetCell(num_pu[m],num_r,num_ru);
								    for (o=0;o<this->DimensionalidadEstados;o++)
								      {
									 pot_aux[o]=0;
									 for (n=0;n<num_r;n++)
									   pot_aux[o]+=(Vecinos->GetComponent(num_qu[num_ru[n]],j+o*this->DimensionalidadEstados)-Vecinos->GetComponent(num_qu[num_ru[n]],(j+1)%this->NumEntidades)+o*this->DimensionalidadEstados)*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
								      }
								    pot+=this->Potencial(pot_aux,p-1+this->CatCliques->GetNumberOfCells()-1);
								 }
							       this->EnergiaNormalizada[u]->SetScalarComponentFromFloat(k,s,t,j*this->L+p-1+this->CatCliques->GetNumberOfCells()-1,pot);
							    }
						       }
						  }
						for (l=0;l<this->V;l++)
						  {
						     pot=-this->AlfaPot->GetValue(l+a)*this->LE[j][u]->GetScalarComponentAsFloat(k,s,t,l);
						     this->EnergiaNormalizada[u]->SetScalarComponentFromFloat(k,s,t,a+l+j*this->L,pot);
						  }
					     }
					   for (k=0;k<this->K;k++)
					     {
						aux=0;
						for (i=0;i<this->L;i++)
						  aux+=this->EnergiaNormalizada[u]->GetScalarComponentAsFloat(k,s,t,i+j*this->L)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],i+j*this->L);
						if (k==0)
						  maxE=aux;
						this->Potencia[u]->SetScalarComponentFromFloat(k,s,t,j,aux);
						if (maxE<aux)
						  maxE=aux;
						//Una posibilidad seria restar estos valores max a los valores de potencia...
					     }
					}
				   }
			      }
			 }
		    }
	       }
	  }
     }
}

vtkStructuredPoints *vtkOptimizaContorno::DevuelveRho()
{
   //   this->PuntObservacion->InsertNextPoint(this->P,this->Q,0);
   this->Monitoriza();
   return this->Rho;
}


void vtkOptimizaContorno::Ejecuta()
{  
   int I,O,U;
   
   for (U=this->Ord;U<2*this->Ord+1;U++)
     {
	for(O=this->Ord;O<2*this->Ord+1;O++)
	  {
	     for (I=this->Ord;I<2*this->Ord+1;I++)
	       {
		  this->PuntObservacion->InsertNextPoint(I,O,U);
		  this->Monitoriza();
		  this->Subredecilla(I,O,U);
	       }	     
	  }       	  
     }
}

void vtkOptimizaContorno::Subredecilla(int I,int O,int U)
{
   int i,k,o,u,j,l,m,b,inda[3],indCampoMedio[2];
   float maxE[2],r;
   double punt[3];
   double Z;
   
   vtkFloatArray *Vecinos=vtkFloatArray::New();
 
   Vecinos->SetNumberOfComponents(this->NumEntidades*this->DimensionalidadEstados);
   if (this->MFA)
     Vecinos->SetNumberOfTuples(this->Vecindario->GetNumberOfPoints()*this->K);
   else     
     Vecinos->SetNumberOfTuples(this->Vecindario->GetNumberOfPoints());
   
   
   b=this->L-this->V;
   for (inda[2]=0;inda[2]<this->DimZ[2];inda[2]++)
     {
	for (inda[1]=0;inda[1]<this->DimZ[1];inda[1]++)
	  {
	     for (inda[0]=0;inda[0]<this->DimZ[0];inda[0]++)
	       {
		  for (u=U+(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],2);u<(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2]+1,2)+this->Ord;u+=this->Ord+1)
		    {
		       for (o=O+(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],1);o<(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1]+1,inda[2],1)+this->Ord;o+=this->Ord+1)
			 {
			    for (i=I+(int)this->Zonas->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],0);i<(int)this->Zonas->GetScalarComponentAsFloat(inda[0]+1,inda[1],inda[2],0)+this->Ord;i+=this->Ord+1)
			      {
				 for (j=0;j<this->NumEntidades;j++)
				   {
				      for (l=1;l<this->Vecindario->GetNumberOfPoints();l++)
					{
					   this->Vecindario->GetPoint(l,punt);
					   for (m=0;m<this->DimensionalidadEstados;m++)
					     {
						if (this->MFA)
						  {
						     for (k=0;k<this->K;k++)
						       {							  
							  r=this->CampoMedio[j]->GetScalarComponentAsFloat((i+(int)punt[0]-this->Ord+this->DimMalla[0])%this->DimMalla[0],(o+(int)punt[1]-this->Ord+this->DimMalla[1])%this->DimMalla[1],(u+(int)punt[2]-this->Ord+this->DimMalla[2])%this->DimMalla[2],k+m*this->K);
							  Vecinos->SetComponent(l+k*this->Vecindario->GetNumberOfPoints(),j+m*this->DimensionalidadEstados,r);
						       }						     
						  }						
						else
						  {						     						       
						     r=this->Deformaciones->GetComponent((int)this->RhoBordes->GetScalarComponentAsFloat(i+(int)punt[0],o+(int)punt[1],u+(int)punt[2],j),m);
						     Vecinos->SetComponent(l,j+m*this->DimensionalidadEstados,r);
						  }						
					     }
					}
				   }
				 for (j=0;j<this->NumEntidades;j++)
				   {
				      for (m=0;m<this->DimensionalidadEstados;m++)
					{
					   for (l=0;l<this->NumEntidades;l++)
					     {
						if (l!=j)
						  {
						     if (this->MFA)
						       {							  
							  for (k=0;k<this->K;k++)
							    {
							       r=this->CampoMedio[j]->GetScalarComponentAsFloat(i-this->Ord,o-this->Ord,u-this->Ord,k+m*this->K);
							       Vecinos->SetComponent(k*this->Vecindario->GetNumberOfPoints(),l+m*this->DimensionalidadEstados,r);							       
							    }
						       }
						     else
						       {							  
							  r=this->Deformaciones->GetComponent((int)this->RhoBordes->GetScalarComponentAsFloat(i,o,u,l),m);
							  Vecinos->SetComponent(0,l+m*this->DimensionalidadEstados,r);
						       }						     
						  }
					     }
					   r=this->Deformaciones->GetComponent(0,m);
					   Vecinos->SetComponent(0,j+m*this->DimensionalidadEstados,r);
					}
				      
				      this->EnE[j]->SetValue(0,this->Energia(inda,j,0,0,Vecinos));
				      for (l=0;l<this->V;l++)
					this->EnE[j]->SetValue(0,this->EnE[j]->GetValue(0)-this->AlfaPot->GetValue(l+b)*this->LE[j][u-this->Ord]->GetScalarComponentAsFloat(0,i-this->Ord,o-this->Ord,l)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j*this->L+l+b));
				      
				      /*if (this->NumEntidades>1 && j!=this->NumEntidades-1)
				       {
				       if (Vecinos->GetComponent(0,j+1)-Vecinos->GetComponent(0,j)<-this->K/3.0)
				       this->EnE[j]->SetValue(0,this->EnE[j]->GetValue(0)-30.0);
				       }
				       else if (this->NumEntidades>1 && j==this->NumEntidades-1)
				       {
				       if (Vecinos->GetComponent(0,j)-Vecinos->GetComponent(0,0)<-this->K/3.0)
				       this->EnE[j]->SetValue(0,this->EnE[j]->GetValue(0)-30.0);
				       }*/
				      if (this->NumEntidades>1)
					this->EnE[j]->SetValue(0,this->EnE[j]->GetValue(0)+this->Energia(inda,j,0,1,Vecinos));
				      this->EnE[j]->SetValue(0,this->EnE[j]->GetValue(0)*this->Beta);
				      
				      maxE[j]=this->EnE[j]->GetValue(0);
				      indCampoMedio[j]=0;
				      for(k=1;k<this->K;k++)
					{
					   for (m=0;m<this->DimensionalidadEstados;m++)
					     {
						r=this->Deformaciones->GetComponent(k,m);
						Vecinos->SetComponent(0,j+m*this->DimensionalidadEstados,r);
					     }
					   this->EnE[j]->SetValue(k,this->Energia(inda,j,k,0,Vecinos));
					   for (l=0;l<this->V;l++)
					     this->EnE[j]->SetValue(k,this->EnE[j]->GetValue(k)-this->AlfaPot->GetValue(l+b)*this->LE[j][u-this->Ord]->GetScalarComponentAsFloat(k,i-this->Ord,o-this->Ord,l)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],j*this->L+l+b));
					   /*if (this->NumEntidades>1 && j!=this->NumEntidades-1)
					    {
					    if (Vecinos->GetComponent(0,j+1)-Vecinos->GetComponent(0,j)<-this->K/3.0)
					    this->EnE[j]->SetValue(k,this->EnE[j]->GetValue(k)-30.0);
					    }
					    else if (this->NumEntidades>1 && j==this->NumEntidades-1)
					    {
					    if (Vecinos->GetComponent(0,j)-Vecinos->GetComponent(0,0)<-this->K/3.0)
					    this->EnE[j]->SetValue(k,this->EnE[j]->GetValue(k)-30.0);
					    }*/
					   if (this->NumEntidades>1)
					     this->EnE[j]->SetValue(k,this->EnE[j]->GetValue(k)+this->Energia(inda,j,k,1,Vecinos));
					   this->EnE[j]->SetValue(k,this->EnE[j]->GetValue(k)*this->Beta);
					   if (maxE[j]<this->EnE[j]->GetValue(k))
					     {
						indCampoMedio[j]=k;
						maxE[j]=this->EnE[j]->GetValue(k);
					     }					   
					}
				   }
				 for (j=0;j<this->NumEntidades;j++)
				   {
				      Z=0;
				      for (k=0;k<this->K;k++)
					{
					   this->Pg->SetValue(k,exp((double)this->EnE[j]->GetValue(k)-maxE[j]));
					   this->CampoMedio[j]->SetScalarComponentFromFloat(i-this->Ord,o-this->Ord,u-this->Ord,k,this->Pg->GetValue(k));
					   Z=Z+this->Pg->GetValue(k);
					   this->Pg->SetValue(k,Z);
					}
				      for (k=0;k<this->K;k++)					
					this->CampoMedio[j]->SetScalarComponentFromFloat(i-this->Ord,o-this->Ord,u-this->Ord,k,this->CampoMedio[j]->GetScalarComponentAsFloat(i-this->Ord,o-this->Ord,u-this->Ord,k)/this->Pg->GetValue(this->K-1));
				      if (this->MFA)
					this->ActualizaRhoBordes(i,o,u,indCampoMedio[j],j);
				      else
					{					   
					   r=(float)(1.0*rand()/(RAND_MAX+1.0));
					   r=r*Z;
					   for (k=0;k<this->K;k++)
					     {
						if (this->Pg->GetValue(k)>r)
						  {
						     this->ActualizaRhoBordes(i,o,u,k,j);
						     break;
						  }
					     }
					}				      
				   }
			      }
			 }
		    }
	       }
	  }
     }
}

float vtkOptimizaContorno::Energia(int *inda,int FD,int k,int aux, vtkFloatArray *Vec)
{
//float vtkOptimizaContorno::Energia(int z,int j,int k,bool FD,float *vecinos)
   int i,l,a,p,m,n,o,ind[5],num_p,num_q,num_r, *num_pu, *num_qu, *num_ru,ver;
   float energia,probCampo;
   
   float * pot;
   pot = new float[this->DimensionalidadEstados];
   a=this->L-this->V;
   //esto es la energia
   p=0;
   energia=0;
   if (aux==0)
     {
	this->CatCliques->InitTraversal();
	this->CatCliques->GetNextCell(num_p,num_pu);
	for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
	  {
	     for (l=0;l<this->NumClNivel->GetValue(i);l++)
	       {
		  p++;
		  this->CatCliques->GetNextCell(num_p,num_pu);
		  for (m=0;m<num_p;m++)
		    {
		       this->Cliques->GetCell(num_pu[m],num_q,num_qu);
		       this->OrdenCliques->GetCell(num_pu[m],num_r,num_ru);
		       if (this->MFA==1)
			 {
			    if (i==2)
			      {				 
				 for (ind[0]=0;ind[0]<this->K;ind[0]++)
				   {				 
				      for (ind[1]=0;ind[1]<this->K;ind[1]++)
					{
					   for (o=0;o<this->DimensionalidadEstados;o++)
					     {
						pot[o]=0;
						ver=0;
						probCampo=1;
						for (n=0;n<num_r;n++)
						  {
						     if (num_qu[num_ru[n]]==0)						       
						       pot[o]+=k*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
						     else
						       {							  
							  pot[o]+=ind[ver]*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
							  probCampo*=Vec->GetComponent(num_qu[num_ru[n]]+ind[ver]*this->Vecindario->GetNumberOfPoints(),FD+o*this->DimensionalidadEstados);
							  ver++;
						       }
						  }
					     }
					   energia+=probCampo*this->Potencial(pot,p-1)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],this->L*FD+p-1);
					}		
				   }
			      }			 		       
			    else
			      {				 
				 for (ind[0]=0;ind[0]<this->K;ind[0]++)
				   {
				      for (o=0;o<this->DimensionalidadEstados;o++)
					{
					   pot[o]=0;
					   for (n=0;n<num_r;n++)
					     {
						if (num_qu[num_ru[n]]==0)						       
						  pot[o]+=k*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
						else
						  {							  
						     pot[o]+=ind[0]*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
						     probCampo=Vec->GetComponent(num_qu[num_ru[n]]+ind[0]*this->Vecindario->GetNumberOfPoints(),FD+o*this->DimensionalidadEstados);
						  }
					     }					
					   energia+=probCampo*this->Potencial(pot,p-1)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],this->L*FD+p-1);
					}		
				   }
			      }
			 }		       
		       else
			 {
			    for (o=0;o<this->DimensionalidadEstados;o++)
			      {
				 pot[o]=0;
				 for (n=0;n<num_r;n++)
				   pot[o]+=Vec->GetComponent(num_qu[num_ru[n]],FD+o*this->DimensionalidadEstados)*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
			      }
			    energia+=this->Potencial(pot,p-1)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],this->L*FD+p-1);
			 }		       
		    }
	       }
	  }
     }   
   else
     {
	this->CatCliques->InitTraversal();
	this->CatCliques->GetNextCell(num_p,num_pu);
	for (i=1;i<this->NumClNivel->GetNumberOfTuples();i++)
	  {
	     for (l=0;l<this->NumClNivel->GetValue(i);l++)
	       {
		  p++;
		  this->CatCliques->GetNextCell(num_p,num_pu);
		  for (m=0;m<num_p;m++)
		    {
		       this->Cliques->GetCell(num_pu[m],num_q,num_qu);
		       this->OrdenCliques->GetCell(num_pu[m],num_r,num_ru);
		       if (this->MFA==1)
			 {
			    if (i==2)
			      {				 
				 for (ind[0]=0;ind[0]<this->K;ind[0]++)
				   {				 
				      for (ind[1]=0;ind[1]<this->K;ind[1]++)
					{
					   for (ind[2]=0;ind[2]<this->K;ind[2]++)
					     {
						for (ind[3]=0;ind[3]<this->K;ind[3]++)
						  {
						     for (ind[4]=0;ind[4]<this->K;ind[4]++)
						       {							  
							  for (o=0;o<this->DimensionalidadEstados;o++)
							    {
							       pot[o]=0;
							       ver=0;
							       probCampo=1;
							       for (n=0;n<num_r;n++)
								 {
								    if (num_qu[num_ru[n]]==0)
								      {									 
									 pot[o]+=(k-ind[ver])*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);									 
									 probCampo*=Vec->GetComponent(ind[ver]*this->Vecindario->GetNumberOfPoints()+num_qu[num_ru[n]],(1+FD)%this->NumEntidades+o*this->DimensionalidadEstados);
									 ver++;

								      }								    
								    else
								      {							  
									 pot[o]+=(ind[ver]-ind[ver+1])*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
									 probCampo*=Vec->GetComponent(num_qu[num_ru[n]]+ind[ver]*this->Vecindario->GetNumberOfPoints(),FD+o*this->DimensionalidadEstados)*Vec->GetComponent(ind[ver+1]*this->Vecindario->GetNumberOfPoints()+num_qu[num_ru[n]],(1+FD)%this->NumEntidades+o*this->DimensionalidadEstados);
									 ver+=2;
								      }
								 }
							    }
							  energia+=probCampo*this->Potencial(pot,p-1)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],this->L*FD+p-1);
						       }
						  }
					     }					   
					}		
				   }
			      }			 		       
			    else
			      {				 
				 for (ind[0]=0;ind[0]<this->K;ind[0]++)
				   {
				      for (ind[1]=0;ind[1]<this->K;ind[1]++)
					{
					   for (ind[2]=0;ind[2]<this->K;ind[2]++)
					     {
						for (o=0;o<this->DimensionalidadEstados;o++)
						  {
						     pot[o]=0;
						     ver=0;
						     probCampo=1;
						     for (n=0;n<num_r;n++)
						       {
							  if (num_qu[num_ru[n]]==0)
							    {								    
							       pot[o]+=(k-ind[ver])*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
							       probCampo*=Vec->GetComponent(ind[ver]*this->Vecindario->GetNumberOfPoints()+num_qu[num_ru[n]],(1+FD)%this->NumEntidades+o*this->DimensionalidadEstados);
							       ver++;
							    }							  
							  else
							    {							  
							       pot[o]+=(ind[ver]-ind[ver+1])*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
							       probCampo*=Vec->GetComponent(num_qu[num_ru[n]]+ind[ver]*this->Vecindario->GetNumberOfPoints(),FD+o*this->DimensionalidadEstados)*Vec->GetComponent(ind[ver+1]*this->Vecindario->GetNumberOfPoints()+num_qu[num_ru[n]],(1+FD)%this->NumEntidades+o*this->DimensionalidadEstados);
							       ver+=2;
							    }
						       }
						  }						
					     }					
					   energia+=probCampo*this->Potencial(pot,p-1)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],this->L*FD+p-1);
					}		
				   }
			      }
			 }
		       else
			 {			    
			    for (o=0;o<this->DimensionalidadEstados;o++)
			      {
				 pot[o]=0;
				 for (n=0;n<num_r;n++)
				   pot[o]+=(Vec->GetComponent(num_qu[num_ru[n]],FD+o*this->DimensionalidadEstados)-Vec->GetComponent(num_qu[num_ru[n]],(1+FD)%this->NumEntidades+o*this->DimensionalidadEstados))*this->Multiplicadores->GetComponent(num_pu[m]+1+n,0);
			      }
			    energia+=this->Potencial(pot,p-1+this->CatCliques->GetNumberOfCells()-1)*this->Parametros->GetScalarComponentAsFloat(inda[0],inda[1],inda[2],this->L*FD+this->CatCliques->GetNumberOfCells()-1+p-1);
			 }		       
		    }
	       }
	  }
     }
   return energia;
}

void vtkOptimizaContorno::ActualizaRhoBordes(int i,int o,int u,int k,int a)
{
   //Habria que pensar en algun modo de reducir algo mas el numero de lineas de esto, aunque es un lio, 
   //porque las condiciones son muchas...
   int x,y,z,lad[3];
   
   this->RhoBordes->SetScalarComponentFromFloat(i,o,u,a,k);
   if (this->PeriodMalla[2]==1 && (u<2*this->Ord || u>this->DimMalla[2]-1))
     {
	if (u<2*this->Ord)
	  lad[2]=1;
	else
	  lad[2]=-1;
	this->RhoBordes->SetScalarComponentFromFloat(i,o,u+lad[2]*this->DimMalla[2],a,k);
     }
   else if (this->PeriodMalla[2]==0 && (u==this->Ord || u==this->DimMalla[2]-1+this->Ord))
     {
	if (u==this->Ord)
	  lad[2]=0;
	else
	  lad[2]=1;
	for (z=0;z<this->Ord;z++)
	  this->RhoBordes->SetScalarComponentFromFloat(i,o,lad[2]*(u+1)+z,a,k);
     }
   
   if (this->PeriodMalla[1]==1 && (o<2*this->Ord || o>this->DimMalla[1]-1))
     {
	if (o<2*this->Ord)
	  lad[1]=1;
	else
	  lad[1]=-1;
	this->RhoBordes->SetScalarComponentFromFloat(i,o+lad[1]*this->DimMalla[1],u,a,k);
	if (this->PeriodMalla[2]==1 && (u<2*this->Ord || u>this->DimMalla[2]-1))
	    {
	       if (u<2*this->Ord)
		 lad[2]=1;
	       else
		 lad[2]=-1;
	       this->RhoBordes->SetScalarComponentFromFloat(i,o+lad[1]*this->DimMalla[1],u+lad[2]*this->DimMalla[2],a,k);
	    }
	else if (this->PeriodMalla[2]==0 && (u==this->Ord || u==this->DimMalla[2]-1+this->Ord))
	  {
	     if (u==this->Ord)
	       lad[2]=0;
	     else
	       lad[2]=1;
	     for (z=0;z<this->Ord;z++)
	       this->RhoBordes->SetScalarComponentFromFloat(i,o+lad[1]*this->DimMalla[1],lad[2]*(u+1)+z,a,k);
	  }
     }
   else if (this->PeriodMalla[1]==0 && (o==this->Ord || o==this->DimMalla[1]-1+this->Ord))
     {
	if (o==this->Ord)
	  lad[1]=0;
	else
	  lad[1]=1;
	for (y=0;y<this->Ord;y++)
	  this->RhoBordes->SetScalarComponentFromFloat(i,lad[1]*(o+1)+y,u,a,k);
	if (this->PeriodMalla[2]==1 && (u<2*this->Ord || u>this->DimMalla[2]-1))
	    {
	       if (u<2*this->Ord)
		 lad[2]=1;
	       else
		 lad[2]=-1;
	       for (y=0;y<this->Ord;y++)
		 this->RhoBordes->SetScalarComponentFromFloat(i,lad[1]*(o+1)+y,u+lad[2]*this->DimMalla[2]+u,a,k);
	    }
	else if (this->PeriodMalla[2]==0 && (u==this->Ord || u==this->DimMalla[2]-1+this->Ord))
	  {
	     if (u==this->Ord)
	       lad[2]=0;
	     else
	       lad[2]=1;
	     for (z=0;z<this->Ord;z++)
	       {
		  for (y=0;y<this->Ord;y++)
		    this->RhoBordes->SetScalarComponentFromFloat(i,lad[1]*(o+1)+y,lad[2]*(u+1)+z,a,k);
	       }
	  }     
     }
   
   if (this->PeriodMalla[0]==1 && (i<2*this->Ord || i>this->DimMalla[0]-1))
     {
	if (i<2*this->Ord)
	  lad[0]=1;
	else
	  lad[0]=-1;
	this->RhoBordes->SetScalarComponentFromFloat(i+lad[0]*this->DimMalla[0],o,u,a,k);
	if (this->PeriodMalla[1]==1 && (o<2*this->Ord || o>this->DimMalla[1]-1))
	  {
	     if (o<2*this->Ord)
	       lad[1]=1;
	     else
	       lad[1]=-1;
	     this->RhoBordes->SetScalarComponentFromFloat(i+lad[0]*this->DimMalla[0],o+lad[1]*this->DimMalla[1],u,a,k);
	     if (this->PeriodMalla[2]==1 && (u<2*this->Ord || u>this->DimMalla[2]-1))
	       {
		  if (u<2*this->Ord)
		    lad[2]=1;
		  else
		    lad[2]=-1;
		  this->RhoBordes->SetScalarComponentFromFloat(i+lad[0]*this->DimMalla[0],o+lad[1]*this->DimMalla[1],u+lad[2]*this->DimMalla[2],a,k);
	       }
	     else if (this->PeriodMalla[2]==0 && (u==this->Ord || u==this->DimMalla[2]-1+this->Ord))
	       {
		  if (u==this->Ord)
		    lad[2]=0;
		  else
		    lad[2]=1;
		  for (z=0;z<this->Ord;z++)
		    this->RhoBordes->SetScalarComponentFromFloat(i+lad[0]*this->DimMalla[0],o+lad[1]*this->DimMalla[1],lad[2]*(u+1)+z,a,k);
	       }
	  }
	else if (this->PeriodMalla[1]==0 && (o==this->Ord || o==this->DimMalla[1]-1+this->Ord))
	  {
	     if (o==this->Ord)
	       lad[1]=0;
	     else
	       lad[1]=1;
	     for (y=0;y<this->Ord;y++)
	       this->RhoBordes->SetScalarComponentFromFloat(i+lad[0]*this->DimMalla[0],lad[1]*(o+1)+y,u,a,k);
	     if (this->PeriodMalla[2]==1 && (u<2*this->Ord || u>this->DimMalla[2]-1))
	       {
		  if (u<2*this->Ord)
		    lad[2]=1;
		  else
		    lad[2]=-1;
		  for (y=0;y<this->Ord;y++)
		    this->RhoBordes->SetScalarComponentFromFloat(i+lad[0]*this->DimMalla[0],lad[1]*(o+1)+y,u+lad[2]*this->DimMalla[2]+u,a,k);
	       }
	     else if (this->PeriodMalla[2]==0 && (u==this->Ord || u==this->DimMalla[2]-1+this->Ord))
	       {
		  if (u==this->Ord)
		    lad[2]=0;
		  else
		    lad[2]=1;
		  for (z=0;z<this->Ord;z++)
		    {
		       for (y=0;y<this->Ord;y++)
			 this->RhoBordes->SetScalarComponentFromFloat(i+lad[0]*this->DimMalla[0],lad[1]*(o+1)+y,lad[2]*(u+1)+z,a,k);
		    }
	       }     
	  }	
     }
  
   else if (this->PeriodMalla[0]==0 && (i==this->Ord || i==this->DimMalla[0]-1+this->Ord))
     {
	if (i==this->Ord)
	  lad[0]=0;
	else
	  lad[0]=1;
	for (x=0;x<this->Ord;x++)
	  this->RhoBordes->SetScalarComponentFromFloat(lad[0]*(i+1)+x,o,u,a,k);
	if (this->PeriodMalla[1]==1 && (o<2*this->Ord || o>this->DimMalla[1]-1))
	  {
	     if (o<2*this->Ord)
	       lad[1]=1;
	     else
	       lad[1]=-1;
	     for (x=0;x<this->Ord;x++)
	       this->RhoBordes->SetScalarComponentFromFloat(lad[0]*(i+1)+x,o+lad[1]*this->DimMalla[1],u,a,k);
	     if (this->PeriodMalla[2]==1 && (u<2*this->Ord || u>this->DimMalla[2]-1))
	       {
		  if (u<2*this->Ord)
		    lad[2]=1;
		  else
		    lad[2]=-1;
		  for (x=0;x<this->Ord;x++)
		    this->RhoBordes->SetScalarComponentFromFloat(lad[0]*(i+1)+x,o+lad[1]*this->DimMalla[1],u+lad[2]*this->DimMalla[2],a,k);
	       }
	     else if (this->PeriodMalla[2]==0 && (u==this->Ord || u==this->DimMalla[2]-1+this->Ord))
	       {
		  if (u==this->Ord)
		    lad[2]=0;
		  else
		    lad[2]=1;
		  for (z=0;z<this->Ord;z++)
		    {
		       for (x=0;x<this->Ord;x++)
			 this->RhoBordes->SetScalarComponentFromFloat(lad[0]*(i+1)+x,o+lad[1]*this->DimMalla[1],lad[2]*(u+1)+z,a,k);
		    }		  
	       }
	  }
	else if (this->PeriodMalla[1]==0 && (o==this->Ord || o==this->DimMalla[1]-1+this->Ord))
	  {
	     if (o==this->Ord)
	       lad[1]=0;
	     else
	       lad[1]=1;
	     for (y=0;y<this->Ord;y++)
	       {
		  for (x=0;x<this->Ord;x++)
		    this->RhoBordes->SetScalarComponentFromFloat(lad[0]*(i+1)+x,lad[1]*(o+1)+y,u,a,k);
	       }	     
	     if (this->PeriodMalla[2]==1 && (u<2*this->Ord || u>this->DimMalla[2]-1))
	       {
		  if (u<2*this->Ord)
		    lad[2]=1;
		  else
		    lad[2]=-1;
		  for (y=0;y<this->Ord;y++)
		    {
		       for (x=0;x<this->Ord;x++)
			 this->RhoBordes->SetScalarComponentFromFloat(lad[0]*(i+1)+x,lad[1]*(o+1)+y,u+lad[2]*this->DimMalla[2]+u,a,k);
		    }		  
	       }
	     else if (this->PeriodMalla[2]==0 && (u==this->Ord || u==this->DimMalla[2]-1+this->Ord))
	       {
		  if (u==this->Ord)
		    lad[2]=0;
		  else
		    lad[2]=1;
		  for (z=0;z<this->Ord;z++)
		    {
		       for (y=0;y<this->Ord;y++)
			 {
			    for (x=0;x<this->Ord;x++)
			      this->RhoBordes->SetScalarComponentFromFloat(lad[0]*(i+1)+x,lad[1]*(o+1)+y,lad[2]*(u+1)+z,a,k);
			 }		       
		    }
	       }     
	  }	
     }
}

float vtkOptimizaContorno::Potencial(float *a,int i)
{
   int j;
   //float phi=0;
   double phi=0;
   for (j=0;j<this->DimensionalidadEstados;j++)
	   phi+=pow(a[j],2);
   phi=sqrt(phi);
   phi=-this->AlfaPot->GetValue(i)*pow(phi/(this->K-1),this->BetaPot->GetValue(i));
   return phi;
}

void vtkOptimizaContorno::ConvierteRhoBordes()
{
   int j,i,t,u,x,y,z;

   this->RhoBordes->Initialize();
   this->RhoBordes->SetNumberOfScalarComponents(this->NumEntidades);
   this->RhoBordes->SetDimensions(this->DimMalla[0]+2*this->Ord,this->DimMalla[1]+2*this->Ord,this->DimMalla[2]+2*this->Ord);
   
   for (i=0;i<this->NumEntidades;i++)
     {
	for (u=this->Ord;u<this->DimMalla[2]+this->Ord;u++)
	  {
	     for (t=this->Ord;t<this->DimMalla[1]+this->Ord;t++)
	       {
		  for (j=this->Ord;j<this->DimMalla[0]+this->Ord;j++)
		    this->RhoBordes->SetScalarComponentFromFloat(j,t,u,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,t-this->Ord,u-this->Ord,i));
		  for (x=0;x<this->Ord;x++)
		    {
		       this->RhoBordes->SetScalarComponentFromFloat(x,t,u,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),t-this->Ord,u-this->Ord,i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,t,u,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+x*this->PeriodMalla[0],t-this->Ord,u-this->Ord,i));
		    }
	       }
	     for (j=this->Ord;j<this->DimMalla[0]+this->Ord;j++)
	       {
		  for (y=0;y<this->Ord;y++)
		    {
		       this->RhoBordes->SetScalarComponentFromFloat(j,y,u,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),u-this->Ord,i));
		       this->RhoBordes->SetScalarComponentFromFloat(j,this->DimMalla[1]+this->Ord+y,u,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,u-this->Ord,i));
		    }		  
	       }
	     for (y=0;y<this->Ord;y++)
	       {
		  for (x=0;x<this->Ord;x++)
		    {
		       this->RhoBordes->SetScalarComponentFromFloat(x,y,u,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),u-this->Ord,i));
		       this->RhoBordes->SetScalarComponentFromFloat(x,this->DimMalla[1]+this->Ord+y,u,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,u-this->Ord,i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,y,u,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+this->PeriodMalla[0]*x,this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),u-this->Ord,i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,this->DimMalla[1]+this->Ord+y,u,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+this->PeriodMalla[0]*x,(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,u-this->Ord,i));
		    }
	       }
	  }
	for (t=this->Ord;t<this->DimMalla[1]+this->Ord;t++)
	  {
	     for (j=this->Ord;j<this->DimMalla[0]+this->Ord;j++)
	       {
		  for (z=0;z<this->Ord;z++)
		    {
		       this->RhoBordes->SetScalarComponentFromFloat(j,t,z,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,t-this->Ord,this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+z),i));
		       this->RhoBordes->SetScalarComponentFromFloat(j,t,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,t-this->Ord,(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		    }		  
	       }
	     for (z=0;z<this->Ord;z++)
	       {		  
		  for (x=0;x<this->Ord;x++)
		    {		       
		       this->RhoBordes->SetScalarComponentFromFloat(x,t,z,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),t-this->Ord,this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+x),i));
		       this->RhoBordes->SetScalarComponentFromFloat(x,t,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),t-this->Ord,(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,t,z,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+this->PeriodMalla[0]*x,t-this->Ord,this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+z),i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,t,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+this->PeriodMalla[0]*x,t-this->Ord,(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		    }
	       }	     
	  }
	for (j=this->Ord;j<this->DimMalla[0]+this->Ord;j++)
	  {
	     for (z=0;z<this->Ord;z++)
	       {		  
		  for (y=0;y<this->Ord;y++)
		    {
		       this->RhoBordes->SetScalarComponentFromFloat(j,y,z,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+z),i));
		       this->RhoBordes->SetScalarComponentFromFloat(j,y,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		       this->RhoBordes->SetScalarComponentFromFloat(j,this->DimMalla[1]+this->Ord+y,z,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+z),i));
		       this->RhoBordes->SetScalarComponentFromFloat(j,this->DimMalla[1]+this->Ord+y,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat(j-this->Ord,(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		    }
	       }	     
	  }	
	for (z=0;z<this->Ord;z++)
	  {		  
	     for (y=0;y<this->Ord;y++)
	       {
		  for (x=0;x<this->Ord;x++)
		    {
		       this->RhoBordes->SetScalarComponentFromFloat(x,y,z,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+z),i));
		       this->RhoBordes->SetScalarComponentFromFloat(x,y,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		       this->RhoBordes->SetScalarComponentFromFloat(x,this->DimMalla[1]+this->Ord+y,z,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+z),i));
		       this->RhoBordes->SetScalarComponentFromFloat(x,this->DimMalla[1]+this->Ord+y,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat(this->PeriodMalla[0]*(this->DimMalla[0]-this->Ord+x),(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,y,z,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+this->PeriodMalla[0]*x,this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+z),i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,y,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+this->PeriodMalla[0]*x,this->PeriodMalla[1]*(this->DimMalla[1]-this->Ord+y),(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,this->DimMalla[1]+this->Ord+y,z,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+this->PeriodMalla[0]*x,(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,this->PeriodMalla[2]*(this->DimMalla[2]-this->Ord+z),i));
		       this->RhoBordes->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,this->DimMalla[1]+this->Ord+y,this->DimMalla[2]+this->Ord+z,i,this->Rho->GetScalarComponentAsFloat((this->DimMalla[0]-1)*(1-this->PeriodMalla[0])+this->PeriodMalla[0]*x,(this->DimMalla[1]-1)*(1-this->PeriodMalla[1])+this->PeriodMalla[1]*y,(this->DimMalla[2]-1)*(1-this->PeriodMalla[2])+this->PeriodMalla[2]*z,i));
		    }
	       }
	  }		
     }
}

void vtkOptimizaContorno::ConvierteRho()
{
   int j,i,t,u,k,ind;
   float max;
	
   for (i=0;i<this->NumEntidades;i++)
     {
	for (u=0;u<this->DimMalla[2];u++)
	  {
	     for (t=0;t<this->DimMalla[1];t++)
	       {
		  for (j=0;j<this->DimMalla[0];j++)
		    if (this->MFA)
		      {
			 max=this->CampoMedio[i]->GetScalarComponentAsFloat(j,t,u,0);
			 ind=0;
			 for (k=1;k<this->K;k++)
			   {
			      if (this->CampoMedio[i]->GetScalarComponentAsFloat(j,t,u,k)>max)
				{				   
				   max=this->CampoMedio[i]->GetScalarComponentAsFloat(j,t,u,k);
				   ind=k;
				}
			   }			 
			 this->Rho->SetScalarComponentFromFloat(j,t,u,i,ind);
		      }
		  else
		    this->Rho->SetScalarComponentFromFloat(j,t,u,i,this->RhoBordes->GetScalarComponentAsFloat(j+this->Ord,t+this->Ord,u+this->Ord,i));
	       }
	  }
     }
}

void vtkOptimizaContorno::EstableceLE(vtkStructuredPoints *Laux, int i, int j)
{
   this->LE[i][j]=vtkStructuredPoints::New();
   this->LE[i][j]->SetScalarTypeToFloat();
   this->LE[i][j]->SetNumberOfScalarComponents(this->V);
   this->LE[i][j]->SetDimensions(this->K,this->DimMalla[0],this->DimMalla[1]);
   this->LE[i][j]->DeepCopy(Laux);
   if (j==0)
     {	
	this->CampoMedio[i]=vtkStructuredPoints::New();
	this->CampoMedio[i]->SetScalarTypeToFloat();
	this->CampoMedio[i]->SetNumberOfScalarComponents(this->K*this->DimensionalidadEstados);
	this->CampoMedio[i]->SetDimensions(this->DimMalla);
     }   
}

void vtkOptimizaContorno::EstableceRhoMed(vtkStructuredPoints *RhoMed)
{
   int j,i,x;
   
   this->RhoMarMed=vtkStructuredPoints::New();
   this->RhoMarMed->SetNumberOfScalarComponents(this->NumEntidades);
   this->RhoMarMed->SetDimensions(this->DimMalla[0]+2*this->Ord,1,1);
   for (i=0;i<this->NumEntidades;i++)
     {
	for (j=this->Ord;j<this->DimMalla[0]+this->Ord;j++)
	  this->RhoMarMed->SetScalarComponentFromFloat(j,0,0,i,RhoMed->GetScalarComponentAsFloat(j-this->Ord,0,0,i));
	for (x=0;x<this->Ord;x++)
	  {
	     this->RhoMarMed->SetScalarComponentFromFloat(x,0,0,i,this->RhoMarMed->GetScalarComponentAsFloat(this->DimMalla[0]+x,0,0,i));
	     this->RhoMarMed->SetScalarComponentFromFloat(this->DimMalla[0]+this->Ord+x,0,0,i,this->RhoMarMed->GetScalarComponentAsFloat(this->Ord+x,0,0,i));
	  }
     }
}

void vtkOptimizaContorno::EstableceRho(vtkStructuredPoints *Rh)
{
   this->Rho->Initialize();
   this->Rho->SetNumberOfScalarComponents(this->NumEntidades);
   this->Rho->SetDimensions(this->DimMalla);
   this->Rho->DeepCopy(Rh);
}



void vtkOptimizaContorno::Monitoriza()
{
   int k;
   char nombre[200];
   
   vtkPolyDataWriter *Escritor=vtkPolyDataWriter::New();
   vtkPolyData *Datos=vtkPolyData::New();

   if (this->UsimagTool==0)
     {	
	//Caso renal:
	k=sprintf(nombre,"%s.MonOpt.vtk",this->Ruta);
   
	//Caso cardiaco:
	//   k=sprintf(nombre,"%s/vtk/MonOpt.vtk",this->Ruta);
   
	Datos->SetPoints(this->PuntObservacion);
	Escritor->SetFileName(nombre);
	Escritor->SetFileTypeToASCII();
	Escritor->SetInput(Datos);
	Escritor->Write();
     }
   
   Datos->Delete();
   Escritor->Delete();
}

void vtkOptimizaContorno::PrintSelf(ostream& os, vtkIndent indent)
{

}
