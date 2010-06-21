/////////////////////////////////////////////////////////////////////////
//	Programa:	Visualization Toolkit								   //
//	Módulo:	vtkValidacion.cxx								   //
//	Descripción:	Recoge métodos para obtener los cortes con cada imagen del modelo
//			obtenido, con objetivo de compararlos con los resultados médicos       //
//	Fecha:	2005/09/09												   //
//	Lenguaje:	C++													   //
//	Autor:	Lucilio Cordero Grande									   //
//		ETSI Telecomunicacion, Universidad de Valladolid			   //
//		Campus Miguel Delibes, Camino del Cementerio, s/n			   //
//		e-mail: lcorgra@troya.lpi.tel.uva.es						   //
/////////////////////////////////////////////////////////////////////////


#include "vtkValidacion.h"
#include "vtkCutter.h"
#include "vtkPlane.h"
#include "vtkMath.h"
#include "vtkTransform.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataWriter.h"
#include "vtkCellArray.h"

vtkValidacion* vtkValidacion::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkValidacion");
  if(ret)
    {
    return (vtkValidacion*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkValidacion;
}

// Constructor
vtkValidacion::vtkValidacion()
{
   this->ConjPunt=vtkPoints::New();
}

//Destructor.
vtkValidacion::~vtkValidacion()
{

}

vtkPolyData *vtkValidacion::Ejecuta()
{
   double UnImag[2][3],OrigImag[2][3],Normal[2][3],aux[2][3],NormalHom[4],OrigHom[4],Orientacion[4],punt[3];
   //De OrigImag, la primera corresponde a IJK y la segunda a RAS; lo mismo para Normal
   int i,j,k,puntos;
	
   vtkPlane *PlanoImplicito=vtkPlane::New();
   vtkPolyData *Auxiliar[2];
   vtkCutter *Cortador=vtkCutter::New();
   vtkMath *Instrumento=vtkMath::New();
   vtkPolyDataWriter *Escribiente=vtkPolyDataWriter::New();
   vtkPolyData *Resultado=vtkPolyData::New();
   vtkPoints *ConjPuntbis=vtkPoints::New();
   vtkCellArray *Celdas=vtkCellArray::New();

   for (i=0;i<2;i++)
     Auxiliar[i]=vtkPolyData::New();
	
   Cortador->SetInput(this->Plantilla);
   Cortador->SetCutFunction(PlanoImplicito);

   for (i=0;i<34;i++)
     {
	for (j=0;j<3;j++)
	  {
	     OrigImag[0][2-j]=this->OrigenImag->GetComponent(i,j);
	     UnImag[0][2-j]=this->UnitXImag->GetComponent(i,j)/(this->Incremento);
	     UnImag[1][2-j]=this->UnitYImag->GetComponent(i,j)/(this->Incremento);
	  }
	Instrumento->Cross(UnImag[0],UnImag[1],Normal[0]);
	NormalHom[3]=0;
	OrigHom[3]=1;
	for (j=0;j<3;j++)
	  {
	     OrigImag[0][j]=pow(-1.0,j)*(OrigImag[0][j]-this->Origen[j])/this->Incremento;
	     NormalHom[j]=Normal[0][j];
	     OrigHom[j]=OrigImag[0][j];
	  }
	this->MatrizIJKaRAS->MultiplyPoint(OrigHom,OrigHom);
	this->MatrizIJKaRAS->MultiplyPoint(NormalHom,NormalHom);
		
      	for (j=0;j<3;j++)
	  {
	     OrigImag[1][j]=OrigHom[j];
	     Normal[1][j]=NormalHom[j];
	  }
	PlanoImplicito->SetOrigin(OrigImag[1]);
	PlanoImplicito->SetNormal(Normal[1]);
	Auxiliar[0]=Cortador->GetOutput();
	Auxiliar[0]->Update();
	puntos=Auxiliar[0]->GetNumberOfPoints();
	
	for (k=0;k<2;k++)
	  {
	     Orientacion[3]=0;
	     for (j=0;j<3;j++)
	       Orientacion[j]=UnImag[k][j];
	     this->MatrizIJKaRAS->MultiplyPoint(Orientacion,Orientacion);
	     for (j=0;j<3;j++)
	       UnImag[k][j]=Orientacion[j];
	  }
	for (j=0;j<3;j++)
	  OrigImag[0][j]=OrigImag[1][j];
	for (k=0;k<puntos;k++)
	  {
	     Auxiliar[0]->GetPoint(k,aux[1]);
	     //if (i==0)
	     //	ConjPuntbis->InsertNextPoint(auxb);
	     for (j=0;j<3;j++)
	       aux[0][j]=aux[1][j]-OrigImag[0][j];

	     punt[0]=floor((Instrumento->Dot(aux[0],UnImag[0]))/(Instrumento->Norm(UnImag[0],3)*Instrumento->Norm(UnImag[0],3)));
	     punt[1]=ceil((Instrumento->Dot(aux[0],UnImag[1]))/(Instrumento->Norm(UnImag[1],3)*Instrumento->Norm(UnImag[1],3)));
	     punt[2]=i;
	     this->ConjPunt->InsertNextPoint(punt);
	     if (i==0 || i==33)
	       {
		  for (j=0;j<3;j++)
		    aux[0][j]=OrigImag[0][j]+punt[0]*UnImag[0][j]+punt[1]*UnImag[1][j];
		  //ConjPuntbis->InsertNextPoint(auxa);
		  if (k==puntos-1)
		    {
		       for (j=0;j<3;j++)
			 punt[j]=OrigImag[0][j];
		       ConjPuntbis->InsertNextPoint(punt);
		       for(j=0;j<3;j++)
			 punt[j]=OrigImag[0][j]+541*UnImag[0][j];
		       ConjPuntbis->InsertNextPoint(punt);
		       for(j=0;j<3;j++)
			 punt[j]=OrigImag[0][j]+371*UnImag[1][j];
		       ConjPuntbis->InsertNextPoint(punt);
		       for(j=0;j<3;j++)
			 punt[j]=OrigImag[0][j]+541*UnImag[0][j]+371*UnImag[1][j];
		       ConjPuntbis->InsertNextPoint(punt);
		       if (i==0)
			 {
			    Celdas->InsertNextCell(4);
			    Celdas->InsertCellPoint(0);
			    Celdas->InsertCellPoint(1);
			    Celdas->InsertCellPoint(3);
			    Celdas->InsertCellPoint(2);
			 }
		       if (i==33)
			 {
			    Celdas->InsertNextCell(4);
			    Celdas->InsertCellPoint(4);
			    Celdas->InsertCellPoint(5);
			    Celdas->InsertCellPoint(7);
			    Celdas->InsertCellPoint(6);
			 }
		    }	
	       }
	  }
	//ConjPuntbis->InsertNextPoint(i,i,i);
	//		ConjPuntbis->InsertNextPoint(i,i,i);
	//		ConjPuntbis->InsertNextPoint(i,i,i);
     }
   Auxiliar[1]->SetPoints(ConjPuntbis);
   Auxiliar[1]->SetPolys(Celdas);
   Resultado->SetPoints(this->ConjPunt);
   Escribiente->SetFileName("resulvalid");
   Escribiente->SetFileTypeToASCII();
   Escribiente->SetInput(Resultado);
   Escribiente->Write();

   Escribiente->SetFileName("resulvalidprueba");
   Escribiente->SetInput(Auxiliar[1]);
   Escribiente->SetFileTypeToASCII();
   Escribiente->Write();

   return Auxiliar[1];	
}

void vtkValidacion::IntroduceMatriza(vtkMatrix4x4 *Pl)
{
   this->MatrizRASaIJK=vtkMatrix4x4::New();
   this->MatrizIJKaRAS=vtkMatrix4x4::New();
   this->MatrizRASaIJK->DeepCopy(Pl);
   this->MatrizRASaIJK->Invert(this->MatrizRASaIJK,this->MatrizIJKaRAS);
}

void vtkValidacion::IntroducePlantilla(vtkPolyData *Plant)
{
   this->Plantilla=vtkPolyData::New();
   this->Plantilla->DeepCopy(Plant);
}

void vtkValidacion::GeneraValores()
{
   this->OrigenImag=vtkFloatArray::New();
   this->UnitXImag=vtkFloatArray::New();
   this->UnitYImag=vtkFloatArray::New();
		
   this->OrigenImag->SetNumberOfComponents(3);
   this->OrigenImag->SetNumberOfTuples(34);
   this->UnitXImag->SetNumberOfComponents(3);
   this->UnitXImag->SetNumberOfTuples(34);
   this->UnitYImag->SetNumberOfComponents(3);
   this->UnitYImag->SetNumberOfTuples(34);	   

   //Esto ahora se haria leyendolo de un fichero con vtkPolyDataReader...
   this->Origen[0]=-64.88655593;
   //this->Origen[0]=-51.87236930;
   //this->Origen[1]=12.35867991;
   this->Origen[1]=22.21972422;
   //this->Origen[2]=-51.0688;
   this->Origen[2]=-50.97569445;
	
   this->Incremento=0.02402132;

   this->OrigenImag->SetTuple3(0,-50.13955777 , 22.58821013 ,-67.36973148);
   this->OrigenImag->SetTuple3(1,-50.15853385 , 22.60207567 ,-67.36316877);
   this->OrigenImag->SetTuple3(2,-50.18056556 , 22.61338511 ,-67.35647292);
   this->OrigenImag->SetTuple3(3,-50.20560489 , 22.62217522 ,-67.34969828);
   this->OrigenImag->SetTuple3(4,-50.23353302 , 22.62855189 ,-67.34291031);
   this->OrigenImag->SetTuple3(5,-50.26415873 , 22.63269329 ,-67.33618503);
   this->OrigenImag->SetTuple3(6,-50.29722762 , 22.63483760 ,-67.32960288);
   this->OrigenImag->SetTuple3(7,-50.33243572 , 22.63526721 ,-67.32324569);
   this->OrigenImag->SetTuple3(8,-50.36944315 , 22.63429347 ,-67.31719592);
   this->OrigenImag->SetTuple3(9,-50.40788423 , 22.63224896 ,-67.31153241);
   this->OrigenImag->SetTuple3(10,-50.44738084,  22.62947930, -67.30633181);
   this->OrigenImag->SetTuple3(11,-50.48754894,  22.62633334, -67.30166238);
   this->OrigenImag->SetTuple3(12,-50.52800863,  22.62315415, -67.29758054);
   this->OrigenImag->SetTuple3(13,-50.56840053,  22.62025987, -67.29412938);
   this->OrigenImag->SetTuple3(14,-50.60839081,  22.61793029 ,-67.29133400);
   this->OrigenImag->SetTuple3(15,-50.64768112,  22.61640540 ,-67.28920244);
   this->OrigenImag->SetTuple3(16,-50.68601501,  22.61587730 ,-67.28772761);
   this->OrigenImag->SetTuple3(17,-50.72318050,  22.61648527 ,-67.28689053);
   this->OrigenImag->SetTuple3(18,-50.75901334,  22.61831766 ,-67.28666127);
   this->OrigenImag->SetTuple3(19,-50.79338859,  22.62141901 ,-67.28699502);
   this->OrigenImag->SetTuple3(20,-50.82621377,  22.62579798 ,-67.28783630);
   this->OrigenImag->SetTuple3(21,-50.85741537,  22.63143445 ,-67.28912249);
   this->OrigenImag->SetTuple3(22,-50.88693203,  22.63828456 ,-67.29078511);
   this->OrigenImag->SetTuple3(23,-50.91472235,  22.64628177 ,-67.29275428);
   this->OrigenImag->SetTuple3(24,-50.94075572,  22.65534798 ,-67.29495937);
   this->OrigenImag->SetTuple3(25,-50.96499723,  22.66541164 ,-67.29733675);
   this->OrigenImag->SetTuple3(26,-50.98740301,  22.67641614 ,-67.29982883);
   this->OrigenImag->SetTuple3(27,-51.00790962,  22.68833306 ,-67.30238425);
   this->OrigenImag->SetTuple3(28,-51.02642970,  22.70116133 ,-67.30495947);
   this->OrigenImag->SetTuple3(29,-51.04286544,  22.71492253 ,-67.30751293);
   this->OrigenImag->SetTuple3(30,-51.05711461,  22.72966083 ,-67.31000945);
   this->OrigenImag->SetTuple3(31,-51.06906852,  22.74544449 ,-67.31242336);
   this->OrigenImag->SetTuple3(32,-51.07862335,  22.76236305 ,-67.31473552);
   this->OrigenImag->SetTuple3(33,-51.08568393,  22.78052492 ,-67.31692821);
   
   this->UnitXImag->SetTuple3(0,-0.00026151  , 0.00015931  , 0.02613643);
   this->UnitXImag->SetTuple3(1,-0.00024915  , 0.00014458  , 0.02613664);
   this->UnitXImag->SetTuple3(2,-0.00023552  , 0.00013090  , 0.02613684);
   this->UnitXImag->SetTuple3(3,-0.00022065  , 0.00011826  , 0.02613703);
   this->UnitXImag->SetTuple3(4,-0.00020462  , 0.00010664  , 0.02613721);
   this->UnitXImag->SetTuple3(5,-0.00018753  , 0.00009598  , 0.02613738);
   this->UnitXImag->SetTuple3(6,-0.00016952  , 0.00008620  , 0.02613754);
   this->UnitXImag->SetTuple3(7,-0.00015071  , 0.00007721  , 0.02613768);
   this->UnitXImag->SetTuple3(8,-0.00013129  , 0.00006890  , 0.02613781);
   this->UnitXImag->SetTuple3(9,-0.00011142  , 0.00006115  , 0.02613792);
   this->UnitXImag->SetTuple3(10,-0.00009127 ,  0.00005384 ,  0.02613801);
   this->UnitXImag->SetTuple3(11,-0.00007104 ,  0.00004682 ,  0.02613809);
   this->UnitXImag->SetTuple3(12,-0.00005089 ,  0.00003998 ,  0.02613815);
   this->UnitXImag->SetTuple3(13,-0.00003099 ,  0.00003318 ,  0.02613819);
   this->UnitXImag->SetTuple3(14,-0.00001151 ,  0.00002632 ,  0.02613821);
   this->UnitXImag->SetTuple3(15,0.00000740  , 0.00001928  , 0.02613822);
   this->UnitXImag->SetTuple3(16,0.00002563  , 0.00001198  , 0.02613821);
   this->UnitXImag->SetTuple3(17,0.00004306  , 0.00000435  , 0.02613819);
   this->UnitXImag->SetTuple3(18,0.00005959  ,-0.00000366  , 0.02613816);
   this->UnitXImag->SetTuple3(19,0.00007514  ,-0.00001209  , 0.02613812);
   this->UnitXImag->SetTuple3(20,0.00008963  ,-0.00002097  , 0.02613807);
   this->UnitXImag->SetTuple3(21,0.00010300  ,-0.00003032  , 0.02613801);
   this->UnitXImag->SetTuple3(22,0.00011518  ,-0.00004013  , 0.02613794);
   this->UnitXImag->SetTuple3(23,0.00012612  ,-0.00005042  , 0.02613787);
   this->UnitXImag->SetTuple3(24,0.00013576  ,-0.00006116  , 0.02613780);
   this->UnitXImag->SetTuple3(25,0.00014403  ,-0.00007236  , 0.02613773);
   this->UnitXImag->SetTuple3(26,0.00015086  ,-0.00008401  , 0.02613766);
   this->UnitXImag->SetTuple3(27,0.00015619  ,-0.00009610  , 0.02613758);
   this->UnitXImag->SetTuple3(28,0.00015992  ,-0.00010867  , 0.02613751);
   this->UnitXImag->SetTuple3(29,0.00016198  ,-0.00012172  , 0.02613744);
   this->UnitXImag->SetTuple3(30,0.00016227  ,-0.00013529  , 0.02613737);
   this->UnitXImag->SetTuple3(31,0.00016072  ,-0.00014940  , 0.02613731);
   this->UnitXImag->SetTuple3(32,0.00015725  ,-0.00016410  , 0.02613724);
   this->UnitXImag->SetTuple3(33,0.00015181  ,-0.00017942  , 0.02613717);
   
   this->UnitYImag->SetTuple3(0,-0.00079643 , -0.02896668 ,  0.00016859);
   this->UnitYImag->SetTuple3(1,-0.00075577 , -0.02896786 ,  0.00015303);
   this->UnitYImag->SetTuple3(2,-0.00071377 , -0.02896900 ,  0.00013865);
   this->UnitYImag->SetTuple3(3,-0.00067042 , -0.02897009 ,  0.00012542);
   this->UnitYImag->SetTuple3(4,-0.00062574 , -0.02897114 ,  0.00011330);
   this->UnitYImag->SetTuple3(5,-0.00057981 , -0.02897214 ,  0.00010223);
   this->UnitYImag->SetTuple3(6,-0.00053274 , -0.02897308 ,  0.00009210);
   this->UnitYImag->SetTuple3(7,-0.00048464 , -0.02897395 ,  0.00008280);
   this->UnitYImag->SetTuple3(8,-0.00043566 , -0.02897475 ,  0.00007419);
   this->UnitYImag->SetTuple3(9,-0.00038595 , -0.02897547 ,  0.00006615);
   this->UnitYImag->SetTuple3(10,-0.00033566 , -0.02897612,   0.00005851);
   this->UnitYImag->SetTuple3(11,-0.00028496 , -0.02897667,   0.00005113);
   this->UnitYImag->SetTuple3(12,-0.00023401 , -0.02897714,   0.00004387);
   this->UnitYImag->SetTuple3(13,-0.00018294 , -0.02897752,   0.00003657);
   this->UnitYImag->SetTuple3(14,-0.00013188 , -0.02897780,   0.00002912);
   this->UnitYImag->SetTuple3(15,-0.00008093 , -0.02897800,   0.00002140);
   this->UnitYImag->SetTuple3(16,-0.00003018 , -0.02897810,   0.00001331);
   this->UnitYImag->SetTuple3(17,0.00002033  ,-0.02897811 ,  0.00000479);
   this->UnitYImag->SetTuple3(18,0.00007056 , -0.02897803 , -0.00000422);
   this->UnitYImag->SetTuple3(19,0.00012051 , -0.02897787 , -0.00001375);
   this->UnitYImag->SetTuple3(20,0.00017020 , -0.02897761 , -0.00002384);
   this->UnitYImag->SetTuple3(21,0.00021965 , -0.02897727 , -0.00003448);
   this->UnitYImag->SetTuple3(22,0.00026889 , -0.02897684 , -0.00004568);
   this->UnitYImag->SetTuple3(23,0.00031798 , -0.02897632 , -0.00005743);
   this->UnitYImag->SetTuple3(24,0.00036697 , -0.02897571 , -0.00006971);
   this->UnitYImag->SetTuple3(25,0.00041592 , -0.02897502 , -0.00008251);
   this->UnitYImag->SetTuple3(26,0.00046486 , -0.02897423 , -0.00009581);
   this->UnitYImag->SetTuple3(27,0.00051383 , -0.02897336 , -0.00010960);
   this->UnitYImag->SetTuple3(28,0.00056286 , -0.02897239 , -0.00012390);
   this->UnitYImag->SetTuple3(29,0.00061195 , -0.02897132 , -0.00013871);
   this->UnitYImag->SetTuple3(30,0.00066111 , -0.02897017 , -0.00015405);
   this->UnitYImag->SetTuple3(31,0.00071033 , -0.02896891 , -0.00016995);
   this->UnitYImag->SetTuple3(32,0.00075961 , -0.02896756 , -0.00018644);
   this->UnitYImag->SetTuple3(33,0.00080893 , -0.02896611 , -0.00020354);
}

void vtkValidacion::PrintSelf(ostream& os, vtkIndent indent)
{
   vtkObject::PrintSelf(os,indent);
}
