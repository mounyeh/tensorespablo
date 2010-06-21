/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMLXB.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/4 14:28:51 $
  Version:   $Revision: 1.1 
=========================================================================*/
#ifndef _itkMLXB_txx
#define _itkMLXB_txx
#include "itkMLXB.h"
#include "itkImageFileWriter.h"

namespace itk
{

template < class TInputImage1, class TInputImage2 >
MLXB< TInputImage1, TInputImage2 >::MLXB()
{
	this->SetNumberOfRequiredInputs( 2 );
	this->SetNumberOfRequiredOutputs( 1 );
	
	m_NLevels = 5;
	m_NIter   = 5;
	m_BlockSize.Fill(5);
	m_SearchSize.Fill(2);
	m_Sample.Fill(6);
	
	m_Sigma.SetSize( TInputImage1::ImageDimension );
	m_Sigma.Fill( 16.0 );
	
	m_Metric  = 6;
	m_Bins    = 20;
	m_Lambda  = 0.5;
	
	m_GaussianRadius.Fill( 2 );
	
	m_Output = NULL;
}

/** Start registration */
template < class TInputImage1, class TInputImage2 >
void MLXB< TInputImage1, TInputImage2 >::Start( void )
{
	//----------------------------------------------------------------------------------------------------
	// Cast input images to internal image type:
	FixedCastFilterPointer  fixedCaster   = FixedCastFilterType::New();
	MovingCastFilterPointer movingCaster  = MovingCastFilterType::New();
	fixedCaster->SetInput(  this->GetInput1() );
	movingCaster->SetInput( this->GetInput2() );
	
	// Subsample casted inputs:
	PyramidPointer          fixedPyramid  = PyramidType::New();
	PyramidPointer          movingPyramid = PyramidType::New();
	fixedPyramid->SetInput(  fixedCaster->GetOutput()  );
	movingPyramid->SetInput( movingCaster->GetOutput() );
	
	// Set number of levels in both pyramids:
	fixedPyramid->SetNumberOfLevels(  m_NLevels );
	movingPyramid->SetNumberOfLevels( m_NLevels );
	
	// Update downsampling scheme:
	fixedPyramid->Update();
	movingPyramid->Update();
	//----------------------------------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------------------
	// The number of refinement level in the interpoaltion:
	unsigned int RL = 5;
	unsigned int NP = 1;
	for( int k=0; k<RL; k++ )
		NP*=2;
	//----------------------------------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------------------
	// Create block-matching filters pair:
	CCFilterPointer         ccFilterf;
	CCFilterPointer         ccFilteri;
	// Create bayesian regularization filters:
	RegFilterPointer        regFilter1f;
	RegFilterPointer        regFilter2f;
	RegFilterPointer        regFilter1i;
	RegFilterPointer        regFilter2i;
	// Create MAP criterion filter:
	MAPFilterPointer        mapFilterf;
	MAPFilterPointer        mapFilteri;
	// Create transforms:
	TransformPointer        transform1      = TransformType::New();
	TransformPointer        transformf2     = TransformType::New();
	TransformPointer        transformi2     = TransformType::New();
	
	// Set parameters of the transforms
	transform1->SetGridRegion(    movingCaster->GetOutput()->GetOrigin(), 
		                          movingCaster->GetOutput()->GetSpacing(),
		                          movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transform1->SetNumberOfPoints( NP, NP, NP );
	
	typename TransformType::ParametersType par1;
	par1.SetSize( transform1->GetNumberOfParameters() );
	transform1->SetParameters( par1 );
	transform1->SetIdentity();
	
	transformf2->SetGridRegion(   movingCaster->GetOutput()->GetOrigin(), 
		                          movingCaster->GetOutput()->GetSpacing(),
		                          movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transformf2->SetNumberOfPoints( 1, 1, 1 );
	transformi2->SetGridRegion(   movingCaster->GetOutput()->GetOrigin(), 
		                          movingCaster->GetOutput()->GetSpacing(),
		                          movingCaster->GetOutput()->GetLargestPossibleRegion().GetSize() );
	transformi2->SetNumberOfPoints( 1, 1, 1 );
	typename TransformType::ParametersType parf2;
	parf2.SetSize( transformf2->GetNumberOfParameters() );
	typename TransformType::ParametersType pari2;
	pari2.SetSize( transformi2->GetNumberOfParameters() );
	transformf2->SetParameters( parf2 );
	transformf2->SetIdentity();
	transformi2->SetParameters( pari2 );
	transformi2->SetIdentity();

	// Create the internal resampling filter for deformations
	ResampleFilterPointer   resampleFilter;
	InterpolatorPointer     interpolator;
	//----------------------------------------------------------------------------------------------------


	for( unsigned int l=0; l<m_NLevels; l++ ){
		this->InvokeEvent( IterationEvent() );
		//------------------------------------------------------------------------------------------------
		// Deform downsampled moving image
		resampleFilter  = ResampleFilterType::New();
		interpolator    = InterpolatorType::New();
		resampleFilter->SetInterpolator( interpolator );
		resampleFilter->SetSize(          movingPyramid->GetOutput( l )->GetLargestPossibleRegion().GetSize() );
		resampleFilter->SetOutputOrigin(  movingPyramid->GetOutput( l )->GetOrigin()                          );
		resampleFilter->SetOutputSpacing( movingPyramid->GetOutput( l )->GetSpacing()                         );
		resampleFilter->SetTransform(     transform1                                                          );
		resampleFilter->SetInput(         movingPyramid->GetOutput( l )                                       );
		resampleFilter->Update();
		//------------------------------------------------------------------------------------------------

		//------------------------------------------------------------------------------------------------
		// Block-matching for this level of resolution:
		// --------------- The forward block-matching:
		ccFilterf = CCFilterType::New();
		// Set parameters:
		ccFilterf->SetMetric(     m_Metric                       );
		ccFilterf->SetBins(       m_Bins                         );
		ccFilterf->SetSample(     m_Sample                       );
		ccFilterf->SetBlockSize(  m_BlockSize                    );
		ccFilterf->SetSearchSize( m_SearchSize                   );
		ccFilterf->SetInput1(     fixedPyramid->GetOutput( l )   );
		ccFilterf->SetInput2(     resampleFilter->GetOutput( )   );
		ccFilterf->Update();
		// --------------- The inverse block-matching:
		ccFilteri = CCFilterType::New();
		// Set parameters:
		ccFilteri->SetMetric(     m_Metric                       );
		ccFilteri->SetBins(       m_Bins                         );
ÿþ
[4676] 19:20:07:447: [INFO ] RPC min calls 1 RPC max calls 500
[4676] 19:20:07:447: [TRACE] ServiceMain: enter
[4720] 19:20:07:447: [INFO ] SPEventHandlerThread: enter (tid=4720)
[4676] 19:20:07:447: [TRACE] calling RpcServerRegisterAuthInfo
[4676] 19:20:07:447: [TRACE] calling RpcServerListen
[4676] 19:20:07:447: [TRACE] RpcServerListen ret'd 1713
[4676] 19:20:07:447: [TRACE] calling RpcServerRegisterIfEx
[4712] 19:20:07:447: [INFO ] SPEventHandlerThread: enter (tid=4712)
[2268] 19:20:07:510: [TRACE] ClientAttach: enter, pid=x494, user='SYSTEM', machine='PABLO1'
[2268] 19:20:07:510: [INFO ] ClientAttach: LookupAccountSidW: User name SYSTEM Domain name NT AUTHORITY
[2268] 19:20:07:510: [INFO ] ClientAttach(SYSTEM): Auth level = 0x6
[2268] 19:20:07:510: [INFO ] PriList: RequestMediaCall=NULL
[2268] 19:20:07:510: [INFO ] LInitialize: calling NewObject ptLineApp 032713A0, pParams->InitContext 800003ff
[2268] 19:20:07:510: [INFO ] LInitialize: NewObject returned hLineApp 103ee
[2268] 19:20:07:510: [INFO ] LInitialize: initializing ptLineApp->InitContext with 800003ff
[2268] 19:20:07:510: [INFO ] LInitialize: initialized ptLineApp->dwAPIVersion with 20000
[2268] 19:20:07:510: [INFO ] ServerInit: NumProviders=4
[2268] 19:20:07:510: [INFO ] ServerInit: ProviderFilename=unimdm.tsp
[2268] 19:20:07:759: [INFO ] ServerInit: u: Calling TSPI_providerInit
[2268] 19:20:07:931: [INFO ] ServerInit: u init'd, dwNumLines=1, dwNumPhones=0
[2268] 19:20:07:9package ClasesPaginacion;

import java.util.Vector;

public class MMU {
	
	private static TablaDePaginas tabla = new TablaDePaginas();
	public static Memoria MP = new Memoria (16 * 1024);
	public static Memoria Disco = new Memoria (128 * 1024);
	
	
	public static int traducirDireccion (int num_pagina) {
		EntradaDeTabla entrada = (EntradaDeTabla) tabla.tabla.elementAt (num_pagina);
		int num_marco;
		
		if (entrada.valido) {
			System.out.println ("La pagina solicitada esta en MP.");
			num_marco = entrada.num_marco;
		}
		
		else {
			
			int sector = entrada.sector;
			Pagina pagina = Disco.leerPagina (sector / 1024);
			num_marco = MP.primerMarcoLibre();
			
			if (num_marco == -1) {				
				num_marco = elegirMarcoLibre();
				MP.cargarPagina (pagina, num_marco);
				Disco.eliminarPagina (sector / 1024);
				entrada.modificar (num_marco, true);				
			}
			
			MP.cargarPagina (pagina, num_marco);
			
		}

		return num_marco;
	}

	
	private static int elegirMarcoLibre () {
		EntradaDeTabla entrada_victima;
		
		int sector_victima, marco_victima;

		do {
			marco_victima = (int) (Math.random() * tabla.tabla.size());
			entrada_victima = (EntradaDeTabla) tabla.tabla.elementAt (marco_victima);
		} while (!entrada_victima.valido);
		
		System.out.print ("Página y marco elegidos al azar (víctima): ");
		System.out.println ("Página " + entrada_victima.num_pagina + ", marco " + entrada_victima.num_marco);

		Pagina pagina_victima = MP.leerPagina (marco_victima);
		sector_victima = 1024 * Disco.cargarPagina (pagina_victima);
		MP.eliminarPagina (marco_victima);
		entrada_victima.modificar (sector_victima, false);
		
		return marco_victima;
	}
	
	
	public static void mostrarTabla () {
		int indice = 0;
		EntradaDeTabla entrada;
		
		System.out.println ("\n#P\t#M\tVALID\tSECTOR\tSIZE");
		
		while (indice < tabla.tabla.size()) {
			entrada = (EntradaDeTabla) tabla.tabla.elementAt(indice);
			
			System.out.print (entrada.num_pagina + "\t");
			
			if (entrada.num_marco != -1)
				System.out.print (entrada.num_marco + "\t");
			else System.out.print ("-\t");

			System.out.print (entrada.valido + "\t");
			
			if (entrada.sector != -1)
				System.out.print (entrada.sector + "\t");
			else System.out.print ("-\t");
			
			System.out.println (entrada.size);
			
			indice++;
		}
		
		System.out.println ();

	}
	
	
	public static Pagina leerPagina (int numero) {
		Pagina pagina;
		EntradaDeTabla entrada = (EntradaDeTabla) tabla.tabla.elementAt (numero);
		int marco;
		
		if (entrada.valido) {
			marco = entrada.num_marco;		
			pagina = MP.leerPagina (marco);
		}
			
		else {
			marco = entrada.sector / 1024;
			pagina = Disco.leerPagina (marco);
		}
		
		return pagina;		
		
	}
	
	
	public static void cargarPagina (char [] array) {
		int num_marco, tam = array.length;
		
		Pagina pagina = new Pagina (array);
		
		num_marco = MP.primerMarcoLibre();
		
		if (num_marco != -1) {
			MP.cargarPagina (pagina, num_marco);
			tabla.nuevaEntradaMP (num_marco, tam);
		}
		
		else {
			num_marco = Disco.primerMarcoLibre();
			Disco.cargarPagina (pagina, num_marco);
			tabla.nuevaEntradaDisco (1024 * num_marco, tam);
		}
		
	}

}



/* Tabla de páginas */
class TablaDePaginas {
	
	protected Vector tabla;
	
	protected int nuevaEntradaMP (int num_marco, int size) {
		EntradaDeTabla nueva_entrada ;
		nueva_entrada = new EntradaDeTabla (tabla.size(), num_marco, true, -1, size);
		tabla.add (nueva_entrada);
		return 0;
	}
	
	protected int nuevaEntradaDisco (int seccion, int size) {
		EntradaDeTabla nueva_entrada ;
		nueva_entrada = new EntradaDeTabla (tabla.size(), -1, false, seccion, size);
		tabla.add (nueva_entrada);
		return 0;
	}
	
	protected TablaDePaginas () {
		tabla = new Vector (10, 5);
		
	}
	
}


/* Entrada de la tabla */
class EntradaDeTabla {
	protected int num_pagina;
	protected int num_marco;
	protected boolean valido;
	protected int sector;
	protected int size;
	
	
	protected void modificar (int numero, boolean val) {
		if (val) {
			num_marco = numero;
			valido = true;
			sector = -1;
		}
		
		else {
			num_marco = -1;
			valido = false;
			sector = numero;
		}
	}
	
	protected EntradaDeTabla (int n_pagina, int n_marco, boolean val, int sec, int tam) {
		num_pagina = n_pagina;
		num_marco = n_marco;
		valido = val;
		sector = sec;
		size = tam;
		
	}
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  