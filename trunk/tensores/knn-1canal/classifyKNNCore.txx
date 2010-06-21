/*=========================================
 Copyright (c) Ruben Cardenes Almeida. August 2006
 This work is licensed under a 
 Creative Commons Atribution 2.5 License 
 http://creativecommons.org/licenses/by/2.5
 ==========================================*/

#ifndef _classifyKNNCore_txx
#define _classifyKNNCore_txx

#include <iostream>
#include "classifyKNNCore.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "voronoiFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <math.h>

namespace itk {


  /**
   * Constructor
   */
template <class TCh1, class TOutput, unsigned int m_ImageDimension>
  classifyKNNCore<TCh1, TOutput, m_ImageDimension>
  :: classifyKNNCore()
  {
    /**
     * Set the number of inputs and outputs for the filter. 
    **/
    this->SetNumberOfRequiredInputs( 1 );
    this->SetNumberOfRequiredOutputs( 1 );
   
    /**
     * A pointer to the output image is created and conected to the output
     **/
    m_K = 1; 
  }// End Constructor

template <class TCh1, class TOutput, unsigned int m_ImageDimension>
void classifyKNNCore<TCh1, TOutput, m_ImageDimension>::GenerateOutputInformation()
{
  // copy output information from input image
  Superclass::GenerateOutputInformation();

  // use user-specified output information
  if ( this->GetInput() == NULL )
    {
    OutputPointer output = this->GetOutput();
    //output->SetLargestPossibleRegion(  this->GetInput()->GetLargestPossibleRegion() );
    }
    
}


template <class TCh1, class TOutput, unsigned int m_ImageDimension>
 void
  classifyKNNCore<TCh1, TOutput, m_ImageDimension>
  ::SetCh1( const TCh1 * channel1 ) 
  {
    SetNthInput(0, const_cast<TCh1 *>( channel1 ));
   
  }


template <class TCh1, class TOutput, unsigned int m_ImageDimension>
 void
  classifyKNNCore<TCh1, TOutput, m_ImageDimension>
  ::SetDiagram( DiagramImageType * diagram ) 
  {
    //SetNthInput(1, const_cast<DiagramImageType *>( diagram ));
    voro = diagram;
  }

template <class TCh1, class TOutput, unsigned int m_ImageDimension>
 void
  classifyKNNCore<TCh1, TOutput, m_ImageDimension>
::GenerateData()
{

	Ch1Pointer ch1 = dynamic_cast<TCh1 *>(this->ProcessObject::GetInput(0));
    
    typename MinimumMaximumImageCalculator<Ch1Type>::Pointer Ch1Calculator=MinimumMaximumImageCalculator<Ch1Type>::New();
    
    Ch1Calculator->SetImage(ch1);
    Ch1Calculator->ComputeMaximum();
    OutputPointer output = this->GetOutput(0);
	output->SetRegions( ch1->GetLargestPossibleRegion() );
    output->Allocate();

    typedef ImageRegionConstIteratorWithIndex<TCh1> Ch1Iterator;

    //    InputIterator inputIt = InputIterator(input,input->GetRequestedRegion());
    Ch1Iterator ch1It = Ch1Iterator(ch1,ch1->GetRequestedRegion());

    DiagramIndexType featureIndex;
    
    //struct NodeType auxNodo;

    //We suppose classes ordered, with the background as zero class 

    typename DiagramImageType::PixelType featureVector;

    int i,count,clasemax;    
    //    int classCount[m_MaxClass+1];    
    //Le sumo uno, aun a costa de
    //tener un elemento mas, para no tener
    //que restar uno cada vez  

    int *classCount = new int[m_Prototipos.size()+1];
    //int classCount[m_Prototipos.size()+1];
       
    for (ch1It.GoToBegin(); !ch1It.IsAtEnd(); ++ch1It){
       featureIndex[0]=(unsigned int)ch1It.Get();
       featureVector=voro->GetPixel(featureIndex);       
       memset(classCount,0,sizeof(int)*(m_Prototipos.size()+1));
       
       for (i=0;i<m_K;i++) {
	      classCount[m_Prototipos[featureVector[i]].clase]++;
       }

       count=0;
       clasemax=0;
       for (i=0;i<m_K;i++) { 
          if (count < classCount[m_Prototipos[featureVector[i]].clase]) {
			count = classCount[m_Prototipos[featureVector[i]].clase];
	        clasemax = m_Prototipos[featureVector[i]].clase;
	      }
	   }

       output->SetPixel(ch1It.GetIndex(),clasemax);
    }//end for (ch1It.GoToBegin()...

    delete[] classCount;
}//end GenerateData



}//end namespace itk


#endif
