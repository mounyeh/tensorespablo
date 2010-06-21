/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDWTFilter.txx$
  Language:  C++
  Date:      $Date: 2005/04/25 $
  Version:   $Revision: 1.0 $
=========================================================================*/
#ifndef _itkWienerFilter_txx
#define _itkWienerFilter_txx

#include "itkWienerFilter.h"
#include "itkObjectFactory.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkCastImageFilter.h"

namespace itk
{
  
template <class TInputImage,class TOutputImage>
WienerFilter<TInputImage,TOutputImage>
::WienerFilter()
{
	m_iterations  = 1;
	m_lambda = 0;
	m_opcion=1;
	m_gamma=0;

	this->SetNumberOfRequiredInputs( 1 );
	this->SetNumberOfRequiredOutputs( 1 );
	
	OutputImagePointer output = OutputImageType::New();
	this->SetNthOutput( 0, output.GetPointer() );
	
        m_Noise = false;
        m_Bias  = false;
}

/**
   * SETINPUT
   **/
  template <class TInputImage, class TOutputImage  >
  void WienerFilter<TInputImage,TOutputImage >
  ::SetInput(const TInputImage * image ) 
  {
    // Process object is not const-correct so the const casting is required.
    SetNthInput(0, const_cast<TInputImage *>(image) );   
  }
  

template<class TInputImage,class TOutputImage>
void WienerFilter<TInputImage,TOutputImage>
::GenerateData( void )
{
  FILE *file;
  unsigned count9;
  int count4;
  //float neigh[27],
  float minVariance,aveVariance,noiseVariance;
  float sum, SNR, localB, localMS, x_0, x_1, *SNR_table;  
  //int voxel;
  float minLocalVariance=0.0, meanLocalVariance=0.0;
  unsigned int posMinLocalVariance;
  unsigned int num_edges;
  unsigned int num_elementos;
  float localMean[6], localVariance[6];
 
 unsigned int ind_edges[6][17]={{9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26},
			     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,14,15,16,17},
			     {1, 2, 4, 5, 7, 8,10,11,14,16,17,19,20,22,23,25,26},
			     {0, 1, 3, 4, 6, 7, 9,10,12,15,16,18,19,21,22,24,25},
			     {0, 1, 2, 3, 4, 5, 9,10,11,12,14,18,19,20,21,22,23},
			     {3, 4, 5, 6, 7, 8,12,14,15,16,17,21,22,23,24,25,26}};	

  InputImageConstPointer input = this->GetInput();
  OutputImagePointer output = this->GetOutput(0);
 

  if (input->GetImageDimension()==3)
    {
      num_edges=6;
      num_elementos=17;
    }
  else if (input->GetImageDimension()==2)
    {
      //float localMean[4], 
      //float localVariance[4];
      num_edges=4;
      num_elementos=5;

     ind_edges[0][0]=0;
     ind_edges[0][1]=1;
     ind_edges[0][2]=2;
     ind_edges[0][3]=3;
     ind_edges[0][4]=5; 
     
     ind_edges[1][0]=0;
     ind_edges[1][1]=1;
     ind_edges[1][2]=3;
     ind_edges[1][3]=6;
     ind_edges[1][4]=7;

     ind_edges[2][0]=1;
     ind_edges[2][1]=2;
     ind_edges[2][2]=3;
     ind_edges[2][3]=7;
     ind_edges[2][4]=8;

     ind_edges[3][0]=3;
     ind_edges[3][1]=5;
     ind_edges[3][2]=6;
     ind_edges[3][3]=7;
     ind_edges[3][4]=8;

    }

 




 
  output->SetLargestPossibleRegion(input->GetRequestedRegion());
  output->SetRequestedRegion(input->GetRequestedRegion());
  output->SetBufferedRegion(input->GetRequestedRegion());
  output->Allocate();
  output->FillBuffer(0);            
    
  FloatImagePointer mean = FloatImageType::New();
  FloatImagePointer variance = FloatImageType::New();
  FloatImagePointer localS = FloatImageType::New();
 
  

  InputImageRegionType region = input->GetLargestPossibleRegion(); 
  mean->SetRegions(region);
  localS->SetRegions(region);
  variance->SetRegions(region);
 
  
  mean->Allocate();
  variance->Allocate();
  localS->Allocate(); 
  
  // Conversion 
  //typedef itk::CastImageFilter<InputImageType, FloatImageType> CasterImageType;
  //CasterImageType::Pointer caster = CasterImageType::New();
  //caster->SetInput(input);
  //FloatImageType::Pointer copy = caster->GetOutput();
  //caster->Update();
  
  FloatImagePointer copy = FloatImageType::New();
  copy->SetRegions(region);
  copy->Allocate();

  typedef itk::ImageRegionConstIterator<InputImageType> IteratorConstType;
  IteratorConstType inputIt = IteratorConstType(input,input->GetRequestedRegion());

  typedef itk::ImageRegionIterator<FloatImageType> IteratorType;
  IteratorType copyIt = IteratorType(copy,copy->GetRequestedRegion());

  for (copyIt.GoToBegin(), inputIt.GoToBegin();!inputIt.IsAtEnd();++inputIt,++copyIt) {
     copyIt.Set((float)inputIt.Get());
  }
  
  myRadiusType radius;
  radius.Fill(1);  
  NeighborhoodIteratorType nit(radius, copy,copy->GetRequestedRegion()); 

  typedef itk::ImageRegionIterator<FloatImageType> FloatIteratorType;
  typedef itk::ImageRegionIterator<OutputImageType> OutputIteratorType;
  OutputIteratorType outIt = OutputIteratorType(output,output->GetRequestedRegion());
  FloatIteratorType meanIt = FloatIteratorType(mean,mean->GetRequestedRegion());
  FloatIteratorType varianceIt = FloatIteratorType(variance,variance->GetRequestedRegion());
  FloatIteratorType localSIt = FloatIteratorType(localS,localS->GetRequestedRegion());
    
  //Remove bias
  if(m_Bias == true){
    SNR_table=(float *)malloc(sizeof(float)*10001);
    file=fopen("SNR.dat","r");
    fread((char *)SNR_table,1,sizeof(float)*10001,file);
    fclose(file);      

    for (meanIt.GoToBegin(), nit.GoToBegin(),localSIt.GoToBegin(),varianceIt.GoToBegin();!nit.IsAtEnd();++nit,++meanIt,++localSIt,++varianceIt) {
      sum = 0;
      if(m_opcion==0 || m_opcion==2) //wiener Homogeneo, wiener Mix
      {
	for (unsigned int i=0;i<nit.Size();i++)
	  {
	  sum += nit.GetPixel(i);
	  }
	meanIt.Set(sum/nit.Size());
	varianceIt.Set(0);
	for (unsigned int i=0;i<nit.Size();i++) {
	  varianceIt.Set((nit.GetPixel(i) - meanIt.Get())*(nit.GetPixel(i) - meanIt.Get())/(nit.Size()-1));
	}
      }
      if(m_opcion==1 || m_opcion==2) //wiener Edge o wiener Mix
      {
	minLocalVariance=0.0;
	for (unsigned int j=0; j<num_edges; j++)
	  {
	    for(unsigned int i=0; i<num_elementos;i++)
	      {
		sum+=nit.GetPixel(ind_edges[j][i]);
	      }
	    localMean[j]=sum/num_elementos;
	    sum=0;
	    for(unsigned int i=0; i<num_elementos;i++)
	      {
		sum+=(nit.GetPixel(ind_edges[j][i])-localMean[j])*(nit.GetPixel(ind_edges[j][i])-localMean[j]);
	      }
	    localVariance[j]=sum/(num_elementos-1);
	    if (localVariance[j]<minLocalVariance || minLocalVariance==0){
	      minLocalVariance=localVariance[j];
	      posMinLocalVariance=j;
	    }
	  }
	if(m_opcion==1)//wiener edge
	  {
	    meanIt.Set(localMean[posMinLocalVariance]);
	    varianceIt.Set(sqrt(minLocalVariance));
	  }
	else if(m_opcion==2)//wiener mix
	  {
	    if(varianceIt.Get()>minLocalVariance){
	      meanIt.Set(localMean[posMinLocalVariance]);
	      varianceIt.Set(sqrt(minLocalVariance));
	    }
	  }
      }

  
      SNR=(meanIt.Get())/(varianceIt.Get());

      if(SNR<=SNR_table[0]) {
	localSIt.Set(0.0);
      }
      else{
	localMS=(meanIt.Get()*meanIt.Get())+(varianceIt.Get()*varianceIt.Get());
	if(SNR>SNR_table[10000]){
	      x_0=10000/1000;
	      x_1=9999/1000;
	      //Considero pendiente 1
	      //  localB=x_0+(SNR-SNR_table[10000])*(x_1-x_0)/(SNR_table[9900]-SNR_table[10000]);
	      localB=x_0+(SNR-SNR_table[10000]);
	}
	else{
	  for(count4=0; count4<10000; count4++){
	    if(SNR<SNR_table[count4]){
	      x_0=(float)count4/1000; 
	      x_1=(float)(count4-1)/1000;
	      localB=x_0+(SNR-SNR_table[count4])*(x_1-x_0)/(SNR_table[count4-1]-SNR_table[count4]);
	      break;
	    }
	  }
	}
	localSIt.Set(localB*sqrt(localMS/(2+(localB*localB))));
      }
    }

   
    for (meanIt.GoToBegin(),localSIt.GoToBegin(),copyIt.GoToBegin();!meanIt.IsAtEnd();++copyIt,++meanIt,++localSIt) {
      copyIt.Set(copyIt.Get()-meanIt.Get()+localSIt.Get());
    }
  }
  //End remove bias
  

  for (count9=0;count9<m_iterations;count9++) {    
    minVariance=0.0;
    aveVariance=0.0;
   
    for (meanIt.GoToBegin(), nit.GoToBegin(),localSIt.GoToBegin(),varianceIt.GoToBegin();!nit.IsAtEnd();++nit,++meanIt,++localSIt,++varianceIt) {
      sum = 0;
     
      if (m_opcion==0 || m_opcion==2) //WienerHom  WienerMix
	{

	  for (unsigned int i=0;i<nit.Size();i++) {
	    sum += nit.GetPixel(i);
	  }
	
	  meanIt.Set(sum/nit.Size());
	  varianceIt.Set(0);
	  for (unsigned int i=0;i<nit.Size();i++) {
	    varianceIt.Set( sqrt( (nit.GetPixel(i) - meanIt.Get())*(nit.GetPixel(i) - meanIt.Get())/(nit.Size()-1)));
	  }
	}
      if(m_opcion==1 || m_opcion==2) //WienerEdge
	{
	  minLocalVariance=0.0;
	for (unsigned int j=0; j<num_edges; j++)
	  {
	    for(unsigned int i=0; i<num_elementos;i++)
	      {
		sum+=nit.GetPixel(ind_edges[j][i]);
	      }
	    localMean[j]=sum/num_elementos;
	    sum=0;
	    for(unsigned int i=0; i<num_elementos;i++)
	      {
		sum+=(nit.GetPixel(ind_edges[j][i])-localMean[j])*(nit.GetPixel(ind_edges[j][i])-localMean[j]);
	      }
	    localVariance[j]=sum/(num_elementos-1);
	    if (localVariance[j]<minLocalVariance || minLocalVariance==0){
	      minLocalVariance=localVariance[j];
	      posMinLocalVariance=j;
	    }
	  }
	if(m_opcion==1)//wienerEdge
	  {
	    meanIt.Set(localMean[posMinLocalVariance]);
	    varianceIt.Set(sqrt(minLocalVariance));
	  }
	else if(m_opcion==2)  //WienerMix
	  {
	    meanLocalVariance=0.0;
	    for(unsigned int i=0; i<num_edges; i++)
	      {
		meanLocalVariance+=localVariance[i]/num_edges;
	      }
	    if(((meanLocalVariance-minLocalVariance)/meanLocalVariance)>m_gamma)
	      {
		meanIt.Set(localMean[posMinLocalVariance]);
		varianceIt.Set(sqrt(minLocalVariance));
	      }
	  }
       }
      if (minVariance==0.0 || minVariance>varianceIt.Get()) {
	minVariance=varianceIt.Get();
      }
      aveVariance+=varianceIt.Get();
    }
    

    aveVariance/=((float)variance->GetLargestPossibleRegion().GetNumberOfPixels());
    noiseVariance=minVariance*(1-m_lambda)+aveVariance*m_lambda;
    std::cout << "noiseVariance " << noiseVariance << std::endl;    
 
    for (meanIt.GoToBegin(),varianceIt.GoToBegin(),copyIt.GoToBegin(),outIt.GoToBegin();!varianceIt.IsAtEnd();++meanIt,++varianceIt,++copyIt,++outIt) { 
	
      if(m_Noise==true){
	varianceIt.Set(varianceIt.Get()-noiseVariance);
	if(varianceIt.Get()<0){
	       varianceIt.Set(0.0);
	}
      }

      
      outIt.Set(static_cast<OutputPixelType>((copyIt.Get()-meanIt.Get())*(varianceIt.Get())/(varianceIt.Get()+noiseVariance)+(meanIt.Get())));
         
    }
    for (copyIt.GoToBegin(),outIt.GoToBegin();!copyIt.IsAtEnd();++copyIt,++outIt) { 
      copyIt.Set(outIt.Get());
    }
  }
}

template <class TInputImage,class TOutputImage>
void WienerFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Iteraciones =   " << m_iterations   << std::endl;
  os << indent << "Lambda = " << m_lambda << std::endl;
}

} // end namespace itk

#endif

