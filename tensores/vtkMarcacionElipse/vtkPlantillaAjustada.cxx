/////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit
//	Módulo:	vtkPlantillaAjustada.cxx
//	Descripción:	A partir de una superficie que aproxima el contorno 3D genera una
//			serie de parámetros útiles en nuestro algoritmo.
//	Fecha:	2004/09/04
//  Última modificación: 2006/05/22
//	Lenguaje:	C++
//	Autor:	Lucilio Cordero Grande
//		ETSI Telecomunicacion, Universidad de Valladolid
//		Campus Miguel Delibes, Camino del Cementerio, s/n
//		e-mail: lcorgra@atenea.lpi.tel.uva.es
//
//////////////////////////////////////////////////////////////////////


#include "vtkPlantillaAjustada.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyDataWriter.h"
#include "vtkCurvatures.h"
#include "vtkPointData.h"
#include "vtkPolyDataReader.h"
#include "vtkStructuredPointsWriter.h"

//-------------------------------------------------------------------------
vtkPlantillaAjustada* vtkPlantillaAjustada::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkPlantillaAjustada");
  if(ret)
    {
    return (vtkPlantillaAjustada*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkPlantillaAjustada;
}

vtkPlantillaAjustada::vtkPlantillaAjustada()
{
   this->Slth=1.0;
   this->J=50;
   this->K=15;

   this->UsimagTool=0;
   this->drmax=18.75;
   this->salto=2;
   
   this->Rinon=vtkPolyData::New();
   this->FichModelo=NULL;
   this->SetFichModelo("Rinon.vtk");
}

vtkPlantillaAjustada::~vtkPlantillaAjustada()
{
}

void vtkPlantillaAjustada::ObtieneParametros(vtkPolyData *Plantilla)
{
   int in,p_mayor,p_menor,JAUX,i,l,p,j,jaux;
   int desorden[4];
   float auxc,auxd,are,Area;
   float *RadAux, *AngAux, *IncrAre, *ArcAux;
   double pi,theta,mintheta;
   double origen[3],normal[3],unit[2][3];
   double cent[3],cenb[3],cenc[3],cend[3];
   double eje[8][3];
   double radi[5],arct[5];
   double eje_coord[2][4];
   
   this->Instrumento=vtkMath::New();
   this->Centro=vtkFloatArray::New();
   this->Param=vtkStructuredPoints::New();
	
   //Se obtienen términos que caracterizan el plano de la imagen
   pi=this->Instrumento->Pi();
   for (l=0;l<3;l++)
     {
	origen[l]=this->Plano->GetElement(l,3);
	normal[l]=this->Plano->GetElement(l,2);
	unit[0][l]=this->Plano->GetElement(l,0);
	unit[1][l]=this->Plano->GetElement(l,1);
	cent[l]=0;
     }
   
   //Se establecen los límites de una ventana ajustada a la zona a segmentar
   for (i=0;i<3;i++)
     Plantilla->GetPoint(4*i,eje[i]);
   this->limite[0]=eje[1][0];
   this->limite[1]=eje[2][0];
   
   this->limite[2]=eje[0][1];
   for (i=0;i<2;i++)
     {
	if (eje[i+1][1]<this->limite[2])
	  this->limite[2]=eje[i+1][1];
     }
   
   for (i=0;i<3;i++)
     Plantilla->GetPoint(2+4*i,eje[i]);
   this->limite[3]=eje[0][1];
   for (i=0;i<2;i++)
     {
	if (eje[i+1][1]>this->limite[3])
	  this->limite[3]=eje[i+1][1];
     }
   
   for (i=0;i<3;i++)
     Plantilla->GetPoint(3+4*i,eje[i]);
   this->limite[4]=eje[0][2];
   for (i=0;i<2;i++)
     {
	if (eje[i+1][2]<this->limite[4])
	  this->limite[4]=eje[i+1][2];
     }
   
   for (i=0;i<3;i++)
     Plantilla->GetPoint(1+4*i,eje[i]);
   this->limite[5]=eje[0][2];
   for (i=0;i<2;i++)
     {
	if (eje[i+1][2]>this->limite[5])
	  this->limite[5]=eje[i+1][2];
     }		
   //Se obtienen límites en profundidad
   //NOTA: EL FACTOR 3 INDICA LA DISTANCIA EN IMÁGENES DE POSICIONES CONSECUTIVAS DEL 
   //CAMPO DE MARKOV EN PROFUNDIDAD; PARA QUE TUVIERAN EL MISMO PESO LA SUAVIDAD RADIAL 
   //Y EN PROFUNDIDAD, ESTE FACTOR DEBERÍA SER APROXIMADAMENTE IGUAL A LA DISTANCIA ENTRE 
   //RAYOS SUCESIVOS.
   p_menor=(int)floor((origen[0]-this->limite[0])/(this->salto*this->Slth)+8/this->salto+1.5);
   p_mayor=(int)floor((this->limite[1]-origen[0])/(this->salto*this->Slth)+8/this->salto+1.5);
   this->P=p_menor+p_mayor-1;

   //Se reserva memoria
   this->Centro->SetNumberOfComponents(3);
   this->Centro->SetNumberOfTuples(this->P);
   this->Param->SetScalarTypeToFloat();
   this->Param->SetNumberOfScalarComponents(2);
   this->Param->SetDimensions(this->P,this->J,1);

   JAUX=this->J*6;
   RadAux = new float[JAUX];
   IncrAre = new float [JAUX];
   ArcAux= new float [JAUX];
   AngAux=new float [JAUX];
	
   for (jaux=0;jaux<JAUX;jaux++)
     AngAux[jaux]=(2.0*pi*jaux)/JAUX;

   for (p=0;p<this->P;p++)
     {
	//Se obtiene el centro y los ejes para cada imagen.
	for (i=0;i<4;i++)
	  Plantilla->GetPoint(i,eje[i]);
	if (p==0)
	  {
	     
	  }
	else if (p<p_menor)
	  {
	     for (j=0;j<4;j++)
	       {
		  Plantilla->GetPoint(4+j,eje[4+j]);
		  eje[j][0]+=((double)p*this->Slth*this->salto/(origen[0]-this->limite[0]))*(eje[4+j][0]-eje[j][0]);
		  for (i=1;i<3;i++)
		    eje[j][i]+=pow((double)p*this->Slth*this->salto/(origen[0]-this->limite[0]),2)*(eje[4+j][i]-eje[j][i]);
	       }
	  }
	else
	  {
	     for (j=0;j<4;j++)
	       {
		  Plantilla->GetPoint(8+j,eje[4+j]);
		  eje[j][0]+=((double)(p-p_menor+1)*this->Slth*this->salto*(eje[4+j][0]-eje[j][0])/(this->limite[1]-origen[0]));
		  for (i=1;i<3;i++)
		    eje[j][i]+=pow((double)(p-p_menor+1)*this->Slth*this->salto/(this->limite[1]-origen[0]),2)*(eje[4+j][i]-eje[j][i]);
	       }
	  }
	for (i=0;i<3;i++)
	  {
	     cent[i]=0;
	     for (j=0;j<4;j++)
	       cent[i]+=eje[j][i];
	     cent[i]/=4.0;
	  }
	
	//Se convierten valores de ejes al plano de la imagen centrado en el centro
	for (i=0;i<4;i++)
	  {
	     radi[i]=this->Instrumento->Distance2BetweenPoints(cent,eje[i]);
	     for (j=0;j<3;j++)
	       eje[i][j]-=cent[j];
	     for (j=0;j<2;j++)
	       eje_coord[j][i]=this->Instrumento->Dot(eje[i],unit[j]);
	     arct[i]=atan2(eje_coord[1][i],eje_coord[0][i]);
	     if (arct[i]<0)
	       arct[i]+=2*pi;
	     desorden[i]=i;
	  }
	for (i=0;i<3;i++)
	  {
	     for (l=i+1;l<4;l++)
	       {
		  if (arct[l]<arct[i])
		    {
		       mintheta=arct[l];
		       arct[l]=arct[i];
		       arct[i]=mintheta;
		       mintheta=radi[l];
		       radi[l]=radi[i];
		       radi[i]=mintheta;
		       in=desorden[l];
		       desorden[l]=desorden[i];
		       desorden[i]=in;
		    }
	       }
	  }
	arct[4]=arct[0]+2*pi;
	radi[4]=radi[0];

	//Se calculan las coordenadas polares de un gran número de puntos de contorno del 
	//template
	Area=0;
	for (jaux=0;jaux<JAUX;jaux++)
	  {
	     mintheta=3*pi;
	     for (i=0;i<4;i++)
	       {
		  theta=AngAux[jaux]-arct[i];
		  if (theta>2*pi)
		    theta=theta-2*pi;
		  if (theta<0)
		    theta=theta+2*pi;
		  if (theta<mintheta)
		    {
		       mintheta=theta;
		       in=i;
		    }
	       }
	     mintheta=pi*mintheta/(2*fabs(arct[in+1]-arct[in]));
	     RadAux[jaux]=sqrt((radi[in]*radi[in+1])/(radi[in+1]*cos(mintheta)*cos(mintheta)+radi[in]*sin(mintheta)*sin(mintheta)));
	     theta=mintheta;
	     mintheta=0;
	     if (jaux!=0)
	       {
		  for(i=0;i<3;i++)
		    {
		       cenb[i]=cent[i]+RadAux[jaux]*(cos(AngAux[jaux])*unit[0][i]+sin(AngAux[jaux])*unit[1][i]);
		       mintheta=mintheta+pow((cenb[i]-cenc[i]),2);
		       cenc[i]=cenb[i];
		    }
		  mintheta=sqrt(mintheta);
		  ArcAux[jaux]=mintheta;
		  Area=Area+ArcAux[jaux];
	       }	
	     else
	       {
		  for (i=0;i<3;i++)
		    {
		       cend[i]=cent[i]+RadAux[jaux]*(cos(AngAux[jaux])*unit[0][i]+sin(AngAux[jaux])*unit[1][i]);
		       cenc[i]=cend[i];
		    }
	       }
			
	     mintheta=0;
	     if (jaux==JAUX-1)
	       {
		  for (i=0;i<3;i++)
		    mintheta=mintheta+pow((cenb[i]-cend[i]),2);
		  mintheta=sqrt(mintheta);
		  ArcAux[0]=mintheta;
		  Area=Area+ArcAux[0];
	       }		
	  }

	//Se introduce el centro
	if (p<p_menor)
	  this->Centro->SetTuple3(p_menor-1-p,cent[0],cent[1],cent[2]);
	else
	  this->Centro->SetTuple3(p,cent[0],cent[1],cent[2]);
		
	//Se calculan los puntos que utiliza el algoritmo
	Area=Area/this->J;
	j=0;
	auxc=0;
	auxd=0;
	for (jaux=0;jaux<JAUX;jaux++)
	  {
	     auxc=RadAux[jaux];
	     are=ArcAux[jaux];
	     auxd=auxd+are;
	     if (auxd>Area)
	       {
		  auxd=auxd-Area;
		  if (p<p_menor)
		    {
		       this->Param->SetScalarComponentFromFloat(p_menor-1-p,j,0,1,auxc);
		       this->Param->SetScalarComponentFromFloat(p_menor-1-p,j,0,0,AngAux[jaux]);
		    }
		  else
		    {
		       this->Param->SetScalarComponentFromFloat(p,j,0,1,auxc);
		       this->Param->SetScalarComponentFromFloat(p,j,0,0,AngAux[jaux]);
		    }
		  j++;
	       }
	  }
	if (p<p_menor)
	  {
	     this->Param->SetScalarComponentFromFloat(p_menor-1-p,this->J-1,0,1,auxc);
	     this->Param->SetScalarComponentFromFloat(p_menor-1-p,this->J-1,0,0,AngAux[JAUX-1]);
	  }
	else
	  {
	     this->Param->SetScalarComponentFromFloat(p,this->J-1,0,1,auxc);
	     this->Param->SetScalarComponentFromFloat(p,this->J-1,0,0,AngAux[JAUX-1]);
	  }
     }
}

vtkStructuredPoints *vtkPlantillaAjustada::GeneraRhoNulo()
{
   int j,p,*puntero;
   vtkStructuredPoints *Rho=vtkStructuredPoints::New();
	
   Rho->SetScalarTypeToInt();
   Rho->SetDimensions(this->J,this->P,1);
	
   for (p=0;p<this->P;p++)
     {
	for (j=0;j<this->J;j++)
	  {
	     puntero=(int *)(Rho->GetScalarPointer(j,p,0));
	     *puntero=(this->K-1)/2;
	  }
     }
   return Rho;
}

vtkPolyData *vtkPlantillaAjustada::ConstruyeModelo(vtkStructuredPoints *Rho)
{
   int j,p,k,i,auxiliar;
   double su;
   double centro[3],punto[3],unit[2][3],*dr;
   char nombre[200];
   char nombrf[200];
   double *aa;
	
   vtkPolyData *Auxiliar[2];
   vtkPoints *Puntos=vtkPoints::New();
   vtkCellArray *Celdas=vtkCellArray::New();
   vtkPolyDataWriter *Escritor=vtkPolyDataWriter::New();
   vtkStructuredPointsWriter *Escritos=vtkStructuredPointsWriter::New();
   vtkCurvatures *CalCur[2];
   vtkDataArray *CamC[2];
   vtkFloatArray *CamI[2];
   
   for (i=0;i<2;i++)
     {
	Auxiliar[i]=vtkPolyData::New();
	CamI[i]=vtkFloatArray::New();
     }
 
   for (j=0;j<2;j++)
     {	
	su=0;
        for (i=0;i<3;i++)
	  {
	     unit[j][i]=this->Plano->GetElement(i,j);
	     su+=pow(unit[j][i],2);
	  }
     }
   
   dr=new double[this->K];
   aa=new double[2];
   j=sprintf(nombre,"%s.Mod.vtk",this->FichModelo);
   j=sprintf(nombrf,"%s.Rho.vtk",this->FichModelo);
   for (k=0;k<this->K;k++)
     dr[k]=(2.0*k*this->drmax)/(this->K-1)-this->drmax;
   k=0;
	
   for (p=0;p<this->P;p++)
     {
	//if (this->MascaraElipsoide->GetComponent(p,0)==1)
	//{
	for (i=0;i<3;i++)
	  centro[i]=this->Centro->GetComponent(p,i);
	for (j=0;j<this->J;j++)
	  {
	     auxiliar=(int)(Rho->GetScalarComponentAsFloat(j,p,0,0));
	     for (i=0;i<3;i++)
	       {
		  //Recuérdese que unity apunta de arriba a abajo de la imagen.
		  punto[i]=centro[i]+(this->Param->GetScalarComponentAsFloat(p,j,0,1)+dr[auxiliar])*(unit[0][i]*cos(this->Param->GetScalarComponentAsFloat(p,j,0,0))+unit[1][i]*sin(this->Param->GetScalarComponentAsFloat(p,j,0,0)));
	       }
	     Puntos->InsertNextPoint(punto);
	     if (j!=0 && k!=0)
	       {
		  Celdas->InsertNextCell(4);
		  Celdas->InsertCellPoint(j+k*this->J);
		  Celdas->InsertCellPoint(j+(k-1)*this->J);
		  Celdas->InsertCellPoint(j-1+(k-1)*this->J);
		  Celdas->InsertCellPoint(j-1+k*this->J);	  
	       }
	  }
	if (k!=0)
	  {
	     Celdas->InsertNextCell(4);
	     Celdas->InsertCellPoint(this->J-1+k*this->J);
	     Celdas->InsertCellPoint(this->J-1+(k-1)*this->J);
	     Celdas->InsertCellPoint((k-1)*this->J);		
	     Celdas->InsertCellPoint(k*this->J);
	  }
	k++;
	//}
     }
	
   Auxiliar[0]->SetPoints(Puntos);
   Auxiliar[0]->SetPolys(Celdas);
   Auxiliar[0]->Update();

   for (i=0;i<2;i++)
     {
	CalCur[i]=vtkCurvatures::New();
	CalCur[i]->SetCurvatureType(i);
	CalCur[i]->SetInput(Auxiliar[0]);
	Auxiliar[1]=CalCur[i]->GetOutput();
	Auxiliar[1]->Update();
	CamC[i]=Auxiliar[1]->GetPointData()->GetScalars()->NewInstance();
	CamC[i]->SetNumberOfComponents(1);
	CamC[i]->SetNumberOfTuples(Auxiliar[1]->GetPointData()->GetNumberOfTuples());
	CamC[i]=Auxiliar[1]->GetPointData()->GetScalars();
	CamI[i]->SetNumberOfComponents(1);
	CamI[i]->SetNumberOfTuples(CamC[i]->GetNumberOfTuples());
     }

   
   CamC[0]->SetName("CurvGauss");
   CamC[1]->SetName("CurvMedia");
   CamI[0]->SetName("IndForma");
   CamI[1]->SetName("IndCurva");
   
   for (i=0;i<CamC[0]->GetNumberOfTuples();i++)
     {
	for (p=0;p<2;p++)
	  aa[p]=CamC[p]->GetTuple1(i);
	CamI[1]->SetTuple1(i,2.0*aa[1]-aa[0]);
	CamI[0]->SetTuple1(i,2.0*atan(-aa[1]/sqrt(pow(aa[1],2.0)-aa[0]))/this->Instrumento->Pi());
	if (CamI[0]->GetTuple1(i) != CamI[0]->GetTuple1(i))
	  CamI[0]->SetTuple1(i,0);
     }
   this->Rinon->Initialize();
   this->Rinon->DeepCopy(Auxiliar[0]);
   for (i=0;i<2;i++)
     {
	p=this->Rinon->GetPointData()->AddArray(CamC[i]);
	p=this->Rinon->GetPointData()->AddArray(CamI[i]);
	aa=CamC[i]->GetRange();
	this->Min[i]=aa[0];
	this->Max[i]=aa[1];
	aa=CamI[i]->GetRange();
	this->Min[i+2]=aa[0];
	this->Max[i+2]=aa[1];
     }
   
   this->Rinon->Update();
   
   if (this->UsimagTool==0)
     {	
	Escritor->SetFileName(nombre);
	Escritor->SetFileTypeToASCII();
	Escritor->SetInput(this->Rinon);
	Escritor->Write();
	Escritos->SetFileName(nombrf);
	Escritos->SetFileTypeToASCII();
	Escritos->SetInput(Rho);
	Escritos->Write();
     }   
   return this->Rinon;
}

vtkPolyData *vtkPlantillaAjustada::LeeModelo()
{
   int j;
   char nombre[200];
   double *rango;
   
   vtkPolyDataReader *Lector=vtkPolyDataReader::New();
   vtkDataArray *CamP[4];

   rango=new double[2];
   j=sprintf(nombre,"%s.Mod.vtk",this->FichModelo);
   Lector->SetFileName(nombre);
   Lector->Update();
   this->Rinon=Lector->GetOutput();
   this->Rinon->Update();
   
   this->Rinon->GetPointData()->SetActiveScalars("CurvGauss");
   CamP[0]=this->Rinon->GetPointData()->GetScalars()->NewInstance();
   CamP[0]->SetNumberOfComponents(1);
   CamP[0]->SetNumberOfTuples(this->Rinon->GetPointData()->GetNumberOfTuples());
   CamP[0]=this->Rinon->GetPointData()->GetScalars();
   this->Rinon->GetPointData()->SetActiveScalars("CurvMedia");
   CamP[1]=this->Rinon->GetPointData()->GetScalars()->NewInstance();
   CamP[1]->SetNumberOfComponents(1);
   CamP[1]->SetNumberOfTuples(this->Rinon->GetPointData()->GetNumberOfTuples());
   CamP[1]=this->Rinon->GetPointData()->GetScalars();
   this->Rinon->GetPointData()->SetActiveScalars("IndForma");
   CamP[2]=this->Rinon->GetPointData()->GetScalars()->NewInstance();
   CamP[2]->SetNumberOfComponents(1);
   CamP[2]->SetNumberOfTuples(this->Rinon->GetPointData()->GetNumberOfTuples());
   CamP[2]=this->Rinon->GetPointData()->GetScalars();
   this->Rinon->GetPointData()->SetActiveScalars("IndCurva");
   CamP[3]=this->Rinon->GetPointData()->GetScalars()->NewInstance();
   CamP[3]->SetNumberOfComponents(1);
   CamP[3]->SetNumberOfTuples(this->Rinon->GetPointData()->GetNumberOfTuples());
   CamP[3]=this->Rinon->GetPointData()->GetScalars();
   for (j=0;j<4;j++)
     {
	rango=CamP[j]->GetRange();
	this->Min[j]=rango[0];
	this->Max[j]=rango[1];
     }
   
   return this->Rinon;
}

void vtkPlantillaAjustada::IntroducePlano(vtkMatrix4x4 *Pl)
{
   this->Plano=vtkMatrix4x4::New();
   
   this->Plano->DeepCopy(Pl);
}

vtkStructuredPoints *vtkPlantillaAjustada::ObtieneParam()
{
   return this->Param;
}

vtkFloatArray *vtkPlantillaAjustada::ObtieneCentro()
{
   return this->Centro;
}

void vtkPlantillaAjustada::EstableceMascaraElipsoide(vtkIntArray *MascElip)
{
   this->MascaraElipsoide=vtkIntArray::New();
   this->MascaraElipsoide->SetNumberOfComponents(1);
   this->MascaraElipsoide->SetNumberOfTuples(this->P);

   this->MascaraElipsoide->DeepCopy(MascElip);
}

void vtkPlantillaAjustada::PrintSelf(ostream& os, vtkIndent indent)
{
}
