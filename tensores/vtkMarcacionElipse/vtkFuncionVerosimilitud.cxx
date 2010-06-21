/////////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit								   //
//	Módulo:	vtkFuncionVerosimilitud.cxx								   //
//	Descripción:	Recoge métodos para estimar los valores de las dos //
//			funciones de verosimilitud empleadas en el algoritmo       //
//	Fecha:	2004/09/09
//  Última modificación: 2006/05/22												   //
//	Lenguaje:	C++													   //
//	Autor:	Lucilio Cordero Grande									   //
//		ETSI Telecomunicacion, Universidad de Valladolid			   //
//		Campus Miguel Delibes, Camino del Cementerio, s/n			   //
//		e-mail: lcorgra@atenea.lpi.tel.uva.es						   //
/////////////////////////////////////////////////////////////////////////


#include "vtkFuncionVerosimilitud.h"
#include "vtkImageResample.h"
#include "vtkImageEuclideanToPolar.h"
#include "vtkImageIdealLowPass.h"
#include "vtkImageFFT.h"
#include "vtkImageRFFT.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkImageGradientMagnitude.h"
#include "vtkImageMedian3D.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataWriter.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkStructuredPointsWriter.h"
#include "vtkPolyDataReader.h"

vtkFuncionVerosimilitud* vtkFuncionVerosimilitud::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkFuncionVerosimilitud");
  if(ret)
    {
    return (vtkFuncionVerosimilitud*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkFuncionVerosimilitud;
}

// Constructor
vtkFuncionVerosimilitud::vtkFuncionVerosimilitud()
{
   int i,j;
   this->Param=vtkStructuredPoints::New();
   this->Centro=vtkFloatArray::New();
   this->UsimagTool=0;
   for (i=0;i<2;i++)
     {
	this->ImEnt[i]=vtkImageData::New();
	this->ImEnt[i]->SetScalarTypeToFloat();
     }
   this->Sub=vtkImageData::New();
   this->MascaraElipsoide=vtkIntArray::New();
   this->Sub->SetScalarTypeToFloat();
   this->K=15;
   this->J=50;
   this->Ng=7;
   //EN REALIDAD DEBE MULTIPLICARSE DRMAX POR EL VALOR DEL MAYOR SEMIEJE
   this->drmax=18.75;
   this->salto=1;
   this->Ruta=NULL;
   this->SetRuta("LR.vtk");
   this->MatrizRASaIJK=vtkMatrix4x4::New();
   this->MatrizIJKaRAS=vtkMatrix4x4::New();
   for (i=0;i<4;i++)
     {
	for (j=0;j<4;j++)
	  this->MatrizRASaIJK->SetElement(i,j,0);
	this->MatrizRASaIJK->SetElement(i,i,1);
     }
   this->MatrizRASaIJK->Invert(this->MatrizRASaIJK,this->MatrizIJKaRAS);
}

//Destructor.
vtkFuncionVerosimilitud::~vtkFuncionVerosimilitud()
{

}

void vtkFuncionVerosimilitud::Ejecuta()
{  
   int p,e;
   char nombre[200];
   
   vtkStructuredPointsWriter *Escritor=vtkStructuredPointsWriter::New();
   //vtkPoints *Puntos=vtkPoints::New();
   //vtkPolyData *Escribeme=vtkPolyData::New();
   //vtkPolyDataWriter *ComprP=vtkPolyDataWriter::New();
	
   this->LR=vtkStructuredPoints::New();
   this->mod=vtkStructuredPoints::New();
   this->Beta=vtkStructuredPoints::New();
   this->Gamma=vtkStructuredPoints::New();
   this->Mask=vtkStructuredPoints::New();
	
   this->MascaraElipsoide->SetNumberOfComponents(1);
   this->MascaraElipsoide->SetNumberOfTuples(this->P);

   this->LR->SetScalarTypeToFloat();
   this->LR->SetNumberOfScalarComponents(2);
   this->LR->SetDimensions(this->K,this->J,this->P);
 	
   this->Beta->SetScalarTypeToFloat();
   this->Gamma->SetScalarTypeToFloat();
   this->Beta->SetNumberOfScalarComponents(200);
   this->Gamma->SetNumberOfScalarComponents(200*15);
   this->Beta->SetDimensions(this->J,this->K,2);
   this->Gamma->SetDimensions(this->J,this->K,2);

   this->mod->SetScalarTypeToInt();
   this->mod->SetNumberOfScalarComponents(1);
   this->mod->SetDimensions(this->J,this->K,2);
 	
   this->Mask->SetScalarTypeToInt();
   this->Mask->SetDimensions(this->M,this->N,1);

   this->dr=new float[this->K];
   this->modmed=new int[this->J];

   this->mu=this->drmax/(this->K-1);

   for (int k=0;k<this->K;k++)
     this->dr[k]=(2.0*k*this->drmax)/(this->K-1)-this->drmax;

   this->CalculaIf();   
   this->CalculaG();
   for (p=0;p<this->P;p++)
     {		
	e=this->CalculaBetaGamma(p);
	if (e==-1)
	  this->MascaraElipsoide->SetComponent(p,0,0);
	else
	  this->MascaraElipsoide->SetComponent(p,0,1);
	this->CalculaLB(p,e);
	this->CalculaLI(p,e);
     }
   if (this->UsimagTool==0)
     {	
	e=sprintf(nombre,"%s.LR.vtk",this->Ruta);
	Escritor->SetFileName(nombre);
	Escritor->SetInput(this->LR);
	Escritor->SetFileTypeToASCII();
        Escritor->Write();
     }   
}

void vtkFuncionVerosimilitud::PrintSelf(ostream& os, vtkIndent indent)
{
   vtkObject::PrintSelf(os,indent);
}

void vtkFuncionVerosimilitud::CalculaIf()
{
   int j,kaux,p,KAUX,k,i,c[2],e,aux[3];
   float x[2],*draux,in[16],LUC[3],RAS[4],IJK[4];
   
   vtkImageData *Ir=vtkImageData::New();
   vtkImageIdealLowPass *LPF=vtkImageIdealLowPass::New();
   vtkImageFFT *TranFour=vtkImageFFT::New();
   vtkImageRFFT *ITranFour=vtkImageRFFT::New();

   this->If=vtkImageData::New();
	
   Ir->SetScalarTypeToFloat();	
   this->If->SetScalarTypeToFloat();	
   this->ImEnt[0]->GetDimensions(aux);

//   Escritor->SetFileName("Moveros.vtk");
//   Escritor->SetFileTypeToAscii();
//   Escritor->SetInput(Escribeme);
   
   KAUX=10*this->K;
   draux=(float *)malloc(KAUX*sizeof(float));
   for (k=0;k<this->K;k++)
     {
	for (kaux=0;kaux<10;kaux++)
	  draux[10*k+kaux]=this->dr[k]+(kaux*2.0*this->mu)/10-this->mu;
     }
	
   Ir->SetDimensions(this->J,KAUX,this->P);
   Ir->SetScalarTypeToFloat();
   for (p=0;p<this->P;p++)
     {
	RAS[3]=1;
	for (j=0;j<this->J;j++)
	  {
	     for (kaux=0;kaux<KAUX;kaux++)
	       {
		  for (i=0;i<3;i++)
		    {
		       //Recuérdese que unity apunta de arriba a abajo de la imagen.
		       RAS[i]=this->Centro->GetComponent(p,i)+(this->Param->GetScalarComponentAsFloat(p,j,0,1)+draux[kaux])*(this->Plano->GetElement(i,0)*cos(this->Param->GetScalarComponentAsFloat(p,j,0,0))+this->Plano->GetElement(i,1)*sin(this->Param->GetScalarComponentAsFloat(p,j,0,0)));
		    }						
		  this->MatrizRASaIJK->MultiplyPoint(RAS,IJK);
		  for (i=0;i<3;i++)
		    {
		       //Recuérdese que unity apunta de arriba a abajo de la imagen.
		       LUC[i]=IJK[i]-this->limites[i];
		    }

		  e=(int)(floor(LUC[2]));
		  for (i=0;i<2;i++)
		    {
		       c[i]=(int)(floor(LUC[i]));
		       x[i]=LUC[i]-c[i];
		       c[i]+=aux[i];
		    }
		  if (e>=1 && e<aux[2]-1)
		    {
		       for (i=0;i<4;i++)
			 {
			    for (k=0;k<4;k++)
			      in[4*i+k]=this->ImEnt[0]->GetScalarComponentAsFloat((c[0]-1+i)%aux[0],(c[1]-1+k)%aux[1],e,0)+this->ImEnt[0]->GetScalarComponentAsFloat((c[0]-1+i)%aux[0],(c[1]-1+k)%aux[1],e-1,0)+this->ImEnt[0]->GetScalarComponentAsFloat((c[0]-1+i)%aux[0],(c[1]-1+k)%aux[1],e+1,0);
			 }
		       in[0]=(1-x[0])*(1-x[1])*in[0]+(1-x[0])*(2-x[1])*in[1]+(1-x[0])*(x[1]+1)*in[2]+(1-x[0])*x[1]*in[3]+(2-x[0])*(1-x[1])*in[4]+(2-x[0])*(2-x[1])*in[5]+(2-x[0])*(x[1]+1)*in[6]+(2-x[0])*x[1]*in[7]+(x[0]+1)*(1-x[1])*in[8]+(x[0]+1)*(2-x[1])*in[9]+(x[0]+1)*(x[1]+1)*in[10]+(x[0]+1)*x[1]*in[11]+x[0]*(1-x[1])*in[12]+x[0]*(2-x[1])*in[13]+x[0]*(x[1]+1)*in[14]+x[0]*x[1]*in[15];
		       in[0]/=48.0;
					
		       /*ina=this->ImEnt->GetScalarComponentAsFloat(c-1,d-1,e,0);
			ina=ina/16.0;*/
		       Ir->SetScalarComponentFromFloat(j,kaux,p,0,(in[0]+.5)/256.0);
		    }
		  else
		    Ir->SetScalarComponentFromFloat(j,kaux,p,0,0);
	       }
	  }
     }
	
   //Es posible que de esta manera que se muestra debajo no consigamos el
   //objetivo de tener un número de puntos por rayo muy superior a K por 
   //medio de interpolación bicúbica.
   TranFour->SetInput(Ir);

   //¿Cuál ha de ser la frecuencia de corte del filtro?
   //Téngase en cuenta que solo filtramos la componente radial.	
   LPF->SetYCutOff(.05);
   //LPF->SetXCutOff(.5);

   //¿CÓMO FILTRAR FORWARD Y BACKWARD? ¿FILTRO DE FASE 0?
   LPF->SetInput(TranFour->GetOutput());
	
   ITranFour->SetInput(LPF->GetOutput());
   ITranFour->Update();
   this->If=ITranFour->GetOutput();
}

void vtkFuncionVerosimilitud::CalculaG()
{
   //AQUÍ HAY QUE METER BAZA PARA CALCULAR BIEN EL GAUSSIANO (QUE HABRÍA QUE 
   //SUBMUESTREAR ACORDE CON LA IMAGEN Y LA IMAGEN MISMA SUBMUESTREADA
   vtkImageGaussianSmooth *Ka=vtkImageGaussianSmooth::New();
   vtkImageGradientMagnitude *Gr=vtkImageGradientMagnitude::New();
   
   this->G=vtkImageData::New();
	
   this->G->SetScalarTypeToFloat();
   
   //Habria que cambiarlo a dimensionalidad 3.
   Ka->SetDimensionality(3);
   Ka->SetRadiusFactors((this->Ng-1)/2.0,(this->Ng-1)/2.0,(this->Ng-1)/2.0);
   Ka->SetStandardDeviation(sqrt((double)this->Ng/2.0),sqrt((double)this->Ng/2.0),sqrt((double)this->Ng/2.0));
   Ka->SetInput(this->Sub);

   //Habria que cambiarlo a dimensionalidad 3.
   Gr->SetDimensionality(3);
   //Gr->SetInput(this->ImEnt[1]);
   Gr->SetInput(Ka->GetOutput());
   Gr->Update();
   this->G=Gr->GetOutput();
}

int vtkFuncionVerosimilitud::CalculaBetaGamma(int p)
{
   int m,n,j,k,i,e,iter,entrada,contador=0,auxil[3],*punteroa;
   float Ifkj,Ifkjm,auxiliar,CEN[3],unit[2][3],RAS[2][4],IJK[2][4];
   double alm,ang,rad,aux,pi,Z[3];
	
   this->g=0;

   for (k=0;k<2;k++)
     {	
	pi=0;
	for (i=0;i<3;i++)
	  unit[k][i]=this->Plano->GetElement(i,k);
     }
   for (j=0;j<this->J;j++)
     {
	this->modmed[j]=0;
	for (k=0;k<this->K;k++)
	  {
	     punteroa=(int *)this->mod->GetScalarPointer(j,k,0);
	     *punteroa=0;
	     punteroa=(int *)this->mod->GetScalarPointer(j,k,1);
	     *punteroa=0;
	  }
     }
   pi=this->Instrumento->Pi();
   for (i=0;i<3;i++)
     {
	//Recuérdese que unity apunta de arriba a abajo de la imagen.
	RAS[0][i]=this->Centro->GetComponent(p,i);
     }
   RAS[0][3]=1;
   this->MatrizRASaIJK->MultiplyPoint(RAS[0],IJK[0]);
   IJK[1][2]=IJK[0][2];
   e=(int)floor((IJK[1][2]-this->limites[2])/2.0);
   this->G->GetDimensions(auxil);
   if (e>=1 && e<auxil[2]-1)
     {
	IJK[1][3]=1;
	for (m=0;m<this->M;m++)
	  {
	     IJK[1][0]=m*2+this->limites[0];
	     for (n=0;n<this->N;n++)
	       {
		  IJK[1][1]=n*2+this->limites[1];
		  this->MatrizIJKaRAS->MultiplyPoint(IJK[1],RAS[1]);
		  //Creo que no hay incompatibilidades entre este ángulo y el de la
		  //obtención de la distribución de ángulos con misma área encerrada. REPASAR.
		  for (i=0;i<3;i++)
		    CEN[i]=RAS[1][i]-RAS[0][i];
		  for (i=0;i<2;i++)
		    Z[i]=this->Instrumento->Dot(CEN,unit[i]);
		  ang=atan2(Z[1],Z[0]);
		  if (ang<0)
		    ang=ang+2*pi;
		  alm=2*pi;
		  iter=0;
		  for (j=0;j<this->J;j++)
		    {
		       aux=fabs(ang-this->Param->GetScalarComponentAsFloat(p,j,0,0));
		       if (aux>pi)
			 aux=2*pi-aux;
		       if (aux<alm)
			 {
			    alm=aux;
			    iter=j;
			 }
		    }
		  rad=this->Instrumento->Norm(CEN);
		  for (k=0;k<this->K;k++)
		    {
		       if(fabs(rad-(this->Param->GetScalarComponentAsFloat(p,iter,0,1)+this->dr[k]))<=this->mu)
			 {
			    auxiliar=this->mod->GetScalarComponentAsFloat(iter,k,1,0);
			    //if (auxiliar<200)
			      //{
				 this->Beta->SetScalarComponentFromFloat(iter,k,0,(int)auxiliar,(float)m);
				 this->Beta->SetScalarComponentFromFloat(iter,k,1,(int)auxiliar,(float)n);
				 auxiliar++;
			      //}
			    punteroa=(int *)this->mod->GetScalarPointer(iter,k,1);
			    *punteroa=(int)auxiliar;
			    this->modmed[iter]=this->modmed[iter]+1;
			    if (rad<=(this->Param->GetScalarComponentAsFloat(p,iter,0,1)+this->dr[k]))
			      {
				 for (i=k;i<this->K;i++)
				   {
				      auxiliar=this->mod->GetScalarComponentAsFloat(iter,i,0,0);
				      //if (auxiliar<3000)
					//{
					   this->Gamma->SetScalarComponentFromFloat(iter,i,0,(int)auxiliar,(float)m);
					   this->Gamma->SetScalarComponentFromFloat(iter,i,1,(int)auxiliar,(float)n);
					   auxiliar++;
					//}
				      punteroa=(int *)this->mod->GetScalarPointer(iter,i,0);
				      *punteroa=(int)auxiliar;
				   }
			      }
			    else
			      {
				 for (i=k+1;i<this->K;i++)
				   {
				      auxiliar=this->mod->GetScalarComponentAsFloat(iter,i,0,0);
				      //if (auxiliar<3000)
					//{
					   this->Gamma->SetScalarComponentFromFloat(iter,i,0,(int)auxiliar,(float)m);
					   this->Gamma->SetScalarComponentFromFloat(iter,i,1,(int)auxiliar,(float)n);
					   auxiliar++;
					//}
				      punteroa=(int *)this->mod->GetScalarPointer(iter,i,0);
				      *punteroa=(int)auxiliar;
				   }
			      }
			    entrada=1;
			    break;
			 }
		       else
			 entrada=0;
		    }
		  
		  if (entrada==1)
		    {
		       Ifkj=pow((double)this->If->GetScalarComponentAsFloat(iter,k*10+5,p,0),2);
		       Ifkj=sqrt(Ifkj+pow((double)this->If->GetScalarComponentAsFloat(iter,k*10+5,p,1),2));
		       if (k==0)
			 Ifkjm=Ifkj;
		       else
			 {
			    //Todo esto, igual que lo anterior, se podría guardar en
			    //sendos vectores, para evitarnos hacer estos cálculos cada
			    //vez.
			    Ifkjm=pow((double)this->If->GetScalarComponentAsFloat(iter,(k-1)*10+5,p,0),2);
			    Ifkjm=sqrt(Ifkjm+pow((double)this->If->GetScalarComponentAsFloat(iter,(k-1)*10+5,p,1),2));
			 }
		       if (Ifkj>Ifkjm)
			 {
			    punteroa=(int *)this->Mask->GetScalarPointer(m,n,0);
			    *punteroa=1;
			    for (i=-(this->salto-1)/2;i<(this->salto+1)/2;i++)
			      this->g+=this->G->GetScalarComponentAsFloat((m+auxil[0])%auxil[0],(n+auxil[1])%auxil[1],e+i,0);
			    contador+=this->salto;
			 }
		       else
			 {
			    punteroa=(int *)this->Mask->GetScalarPointer(m,n,0);
			    *punteroa=0;
			 }
		    }
		  else
		    {
		       punteroa=(int *)this->Mask->GetScalarPointer(m,n,0);
		       *punteroa=0;
		    }
	       }
	  }
	this->g=this->g/contador;
	return e;
     }
   else
     return -1;
}

void vtkFuncionVerosimilitud::CalculaLB(int p,int e)
{
   int a,k,j,m,n,i,l,auxiliar,aux[3];
   float B,auxiliarbis;
   
   this->G->GetDimensions(aux);
   if (e!=-1)
     {
	for (j=0;j<this->J;j++)
	  {
	     for (k=0;k<this->K;k++)
	       {
		  auxiliar=(int)this->mod->GetScalarComponentAsFloat(j,k,1,0);
		  auxiliarbis=0;
		  for (i=0;i<auxiliar;i++)
		    {
		       m=(int)(this->Beta->GetScalarComponentAsFloat(j,k,0,i));
		       n=(int)(this->Beta->GetScalarComponentAsFloat(j,k,1,i));
		       a=(int)(this->Mask->GetScalarComponentAsFloat(m,n,0,0));
		       
		       if (a!=0 && this->g!=0)
			 {
			    B=this->salto;
			    for (l=-(this->salto-1)/2;l<(this->salto+1)/2;l++)
			      B-=exp(-this->G->GetScalarComponentAsFloat((m+aux[0])%aux[0],(n+aux[1])%aux[1],e+l,0)*log(4.0/3.0)/this->g);
			    auxiliarbis-=B/(this->salto*auxiliar);
			 }
		    }
		  this->LR->SetScalarComponentFromFloat(k,j,p,1,auxiliarbis);
	       }
	  }
     }
   else
     {
	for (k=0;k<this->K;k++)
	  {
	     for (j=0;j<this->J;j++)
	       this->LR->SetScalarComponentFromFloat(k,j,p,1,0);
	  }
     }
   
}

void vtkFuncionVerosimilitud::CalculaLI(int p,int e)
{
   int m,n,j,i,h,k,r,l,max,auxiliara,auxilio[3];
   float aux,auxbis,media,varianza,alfa[2],maxim,minim,auxiliar,auxiliarbis, contador;
   
   vtkFloatArray *Aux=vtkFloatArray::New();

   Aux->SetNumberOfComponents(1);
   
   this->Sub->GetDimensions(auxilio);
   max=this->modmed[0];
   Aux->SetNumberOfTuples(this->salto*max);
   if (e!=-1)
     {
	for (j=0;j<this->J;j++)
	  {
	     if (this->modmed[j]>max)
	       {
		  max=this->modmed[j];
		  Aux->SetNumberOfTuples(this->salto*max);
	       }
	     media=0;
	     varianza=0;
	     contador=0;
	     for (k=0;k<this->K;k++)
	       {
		  auxiliara=(int)this->mod->GetScalarComponentAsFloat(j,k,1,0);
		  for (i=0;i<auxiliara;i++)
		    {
		       m=(int)(this->Beta->GetScalarComponentAsFloat(j,k,0,i));
		       n=(int)(this->Beta->GetScalarComponentAsFloat(j,k,1,i));
		       for (l=0;l<this->salto;l++)
			 {
			    auxiliarbis=(this->Sub->GetScalarComponentAsFloat((m+auxilio[0])%auxilio[0],(n+auxilio[1])%auxilio[1],e+l-(this->salto-1)/2,0)+.5)/256.0;
			    if (auxiliarbis<0.95 && auxiliarbis>0.05)
			      {
				 Aux->SetComponent((int)contador,0,auxiliarbis);
				 contador++;
			      }
			 }
		    }
	       }
	     for (i=0;i<contador;i++)
	       {
		  for (h=i;h<contador;h++)
		    {
		       if (Aux->GetComponent(i,0)>Aux->GetComponent(h,0))
			 {
			    aux=Aux->GetComponent(h,0);
			    Aux->SetComponent(h,0,Aux->GetComponent(i,0));
			    Aux->SetComponent(i,0,aux);
			 }
		    }
	       }
	     
		   if (fmod(contador,2)==0)
	       {
		  for (i=0;i<contador/2;i++)
		    {
		       media=media+Aux->GetComponent(i,0);
		       varianza=varianza+pow((double)Aux->GetComponent(i,0),2);
		    }
		  media=2*media/contador;
		  varianza=2*varianza/contador;
	       }
	     else
	       {
		  for (i=0;i<(contador-1)/2;i++)
		    {
		       media=media+Aux->GetComponent(i,0);
		       varianza=varianza+pow((double)Aux->GetComponent(i,0),2);
		    }
		  media=2*media/(contador-1);
		  varianza=2*varianza/(contador-1);
	       }
	     if (contador>9)
	       {
		  varianza=varianza-pow((double)media,2);
		  if (varianza!=0)
		    {
		       aux=-1+(1-media)*media/varianza;
		       alfa[0]=media*aux;
		       alfa[1]=(1-media)*aux;
		    }
		  else
		    {
		       alfa[0]=1;
		       alfa[1]=1;
		    }
	       }
	     else
	       {
		  alfa[0]=1;
		  alfa[1]=1;
	       }
	     maxim=0;
	     minim=0;
	     for (k=this->K-1;k>=0;k--)
	       {
		  auxiliara=(int)this->mod->GetScalarComponentAsFloat(j,k,0,0);
		  auxiliar=0;
		  if (auxiliara>6)
		    {
		       contador=0;
		       for (i=0;i<auxiliara;i++)
			 {
			    m=(int)(this->Gamma->GetScalarComponentAsFloat(j,k,0,i));
			    n=(int)(this->Gamma->GetScalarComponentAsFloat(j,k,1,i));
			    for(r=-(this->salto-1)/2;r<(this->salto+1)/2;r++)
			      {
				 auxiliarbis=(this->Sub->GetScalarComponentAsFloat((m+auxilio[0])%auxilio[0],(n+auxilio[1])%auxilio[1],e+r,0)+.5)/256.0;
				 if (auxiliarbis<0.95 && auxiliarbis>0.05)
				   {
				      contador++;
				      auxiliar+=((1-alfa[1])*log(1-auxiliarbis)+(1-alfa[0])*log(auxiliarbis));
				   }
			      }	
			    
			 }
		       if (contador!=0)
			 auxiliar/=(float)contador;
		       this->LR->SetScalarComponentFromFloat(k,j,p,0,auxiliar);
		       if (k==this->K-1)
			 {
			    maxim=auxiliar;
			    minim=auxiliar;
			 }		       
		       if(auxiliar>maxim)
			 maxim=auxiliar;
		       if(auxiliar<minim)
			 minim=auxiliar;
		    }
		  else if (k!=this->K-1)
		    this->LR->SetScalarComponentFromFloat(k,j,p,0,this->LR->GetScalarComponentAsFloat(k+1,j,p,0));
		  
		  else
		    this->LR->SetScalarComponentFromFloat(k,j,p,0,0);
	       }
	     if (maxim!=minim)
	       {
		  for (k=this->K-1;k>=0;k--)
		    {
		       auxbis=this->LR->GetScalarComponentAsFloat(k,j,p,0);
		       auxbis=(auxbis-maxim)/(maxim-minim);
		       this->LR->SetScalarComponentFromFloat(k,j,p,0,auxbis);
		    }
	       }
	  }
     }
   else
     {
	for (k=0;k<this->K;k++)
	  {
	     for (j=0;j<this->J;j++)
	       this->LR->SetScalarComponentFromFloat(k,j,p,0,0);
	  }
     }
}

void vtkFuncionVerosimilitud::IntroduceParam(vtkStructuredPoints *Par)
{
   this->Param->SetScalarTypeToFloat();
   this->Param->SetNumberOfScalarComponents(2);
   this->Param->SetDimensions(this->P,this->J,1);
   this->Param->DeepCopy(Par);
}

void vtkFuncionVerosimilitud::IntroduceCentro(vtkFloatArray *Cent)
{
   this->Centro->SetNumberOfComponents(3);
   this->Centro->SetNumberOfTuples(this->P);
   this->Centro->DeepCopy(Cent);
}

void vtkFuncionVerosimilitud::IntroduceImagen(vtkImageData *ImagEnt, int FILT)
{
   int m,n,o,aux[3];
   float limit[2][4];
   
   vtkImageResample *Re=vtkImageResample::New();
	
   ImagEnt->GetDimensions(aux);
	
   this->Mbis=aux[0];
   this->Nbis=aux[1];
   for (m=0;m<2;m++)
     {
	limit[m][0]=this->limite[m]+(2*m-1)*9;
	limit[m][1]=this->limite[2+m]+(2*m-1)*this->drmax;
	limit[m][2]=this->limite[4+m]+(2*m-1)*this->drmax;
	limit[m][3]=1;
	this->MatrizRASaIJK->MultiplyPoint(limit[m],limit[m]);
	limit[m][0]+=(1-2*m)*((this->Ng-1)/2.0+10);
	limit[m][1]+=(2*m-1)*((this->Ng-1)/2.0+10);
	limit[m][2]+=(2*m-1)*2;
     }
	
   this->limites[0]=(int)floor(limit[1][0])-20;
   this->limites[3]=(int)ceil(limit[0][0])+20;
   this->limites[1]=(int)floor(limit[0][1])-35;
   this->limites[4]=(int)ceil(limit[1][1])+35;
   this->limites[2]=(int)floor(limit[0][2]);
   this->limites[5]=(int)ceil(limit[1][2]);
	
   //if (this->limites[2]<0)
   this->limites[2]=0;
   
   //if (this->limites[5]>=aux[2])
   this->limites[5]=aux[2]-1;
   this->limites[0]=0;
   this->limites[1]=0;
   this->limites[3]=aux[0]-1;
   this->limites[4]=aux[1]-1;
   
   this->ImEnt[FILT]->SetDimensions(this->limites[3]-this->limites[0]+1,this->limites[4]-this->limites[1]+1,this->limites[5]-this->limites[2]+1);

   for (o=this->limites[2];o<=this->limites[5];o++)
     {
	for (n=this->limites[1];n<=this->limites[4];n++)
	  {
	     for (m=this->limites[0];m<=this->limites[3];m++)
	       this->ImEnt[FILT]->SetScalarComponentFromFloat(m-this->limites[0],n-this->limites[1],o-this->limites[2],0,ImagEnt->GetScalarComponentAsFloat((m+aux[0])%aux[0],(n+aux[1])%aux[1],(o+aux[2])%aux[2],0));
	  }
     }
   if (FILT==0)
     {
	Re->SetDimensionality(3);
	Re->SetAxisMagnificationFactor(0,0.5);
	Re->SetAxisMagnificationFactor(1,0.5);
	Re->SetAxisMagnificationFactor(2,0.5);
	Re->SetInput(this->ImEnt[0]);
	Re->Update();
	this->Sub=Re->GetOutput();
	this->Sub->GetDimensions(aux);
	this->M=aux[0];
	this->N=aux[1];
     }
}


void vtkFuncionVerosimilitud::EscribeImagen()
{
   //Un escritor TIFF o JPEG (para lo cual hay que hacer un ImageCast).
}

vtkStructuredPoints *vtkFuncionVerosimilitud::DevuelveLR()
{
   return this->LR;
}

vtkIntArray *vtkFuncionVerosimilitud::DevuelveMascaraElipsoide()
{
   return this->MascaraElipsoide;
}

void vtkFuncionVerosimilitud::IntroduceMatriza(vtkMatrix4x4 *Pl)
{
//   int e,i;
//   double punt[3];
//   char nombre[200];
//   vtkPolyDataReader *Lector=vtkPolyDataReader::New();
//   vtkPolyDataWriter *Escritor=vtkPolyDataWriter::New();
//   vtkPoints *Puntos=vtkPoints::New();
//   vtkPolyData *Poli=vtkPolyData::New();
//   e=sprintf(nombre,"%s.RASaIJK.vtk",this->Ruta);
   //e=sprintf(nombre,"%s",this->Ruta);
/*   if (this->UsimagTool==1)
     {
	Lector->SetFileName(nombre);
	Poli=Lector->GetOutput();
        Poli->Update();
	for (i=0;i<4;i++)
	  {
	     Poli->GetPoint(i,punt);
	     for (e=0;e<3;e++)
	       this->MatrizRASaIJK->SetElement(i,e,punt[e]);
	  }
	Poli->GetPoint(i,punt);
	for (e=0;e<3;e++)
	  this->MatrizRASaIJK->SetElement(e,3,punt[e]);
	this->MatrizRASaIJK->SetElement(3,3,1);
	
     }*/
//   else
//    {
	this->MatrizRASaIJK->DeepCopy(Pl);
/*	for (i=0;i<4;i++)
	  {	     
	     for (e=0;e<3;e++)
	       punt[e]=this->MatrizRASaIJK->GetElement(i,e);
	     Puntos->InsertNextPoint(punt);
	  }
	for (e=0;e<3;e++)
	  punt[e]=this->MatrizRASaIJK->GetElement(e,3);
	Puntos->InsertNextPoint(punt);
	Escritor->SetFileName(nombre);
	Poli->SetPoints(Puntos);
	Escritor->SetInput(Poli);
	Escritor->Write();     	
    }*/
   this->MatrizRASaIJK->Invert(this->MatrizRASaIJK,this->MatrizIJKaRAS);   
}

void vtkFuncionVerosimilitud::IntroducePlano(vtkMatrix4x4 *Pl)
{
   this->Plano=vtkMatrix4x4::New();
   this->Plano->DeepCopy(Pl);
}

