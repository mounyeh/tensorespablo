/*
 *  itkTractography.txx
 *
 */
// 
#include "itkTractography.h"

#ifndef __itkTractography_txx
#define __itkTractography_txx

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include <itkSymmetricSecondRankTensor.h>

namespace itk {

   /**
   * Constructor
   */
  
  template <class TInput, class TOutput> 
  Tractography<TInput, TOutput> ::Tractography(){ 

	this->SetCurvatureThreshold(90);
	this->SetStepLength(0.1);
	this->SetFaThreshold(0.15);
	this->SetMinimumFiberLength(2);
	this->SetFibersPerVoxel(10);
	this->SetConnection(false);

	m_StreamlineList.clear();
	Superclass::SetNumberOfThreads(1);
/*	this->SetNumberOfRequiredInputs( 1 );
    this->SetNumberOfRequiredOutputs( 1 );
*/
  };
   
   
   
    /**
   * SETINPUT
   **/
  /*template <class TInput, class TOutput  >
  void Tractography<TInput,TOutput >
  ::SetInput(const  TInput * domain ) 
  {
    // Process object is not const-correct so the const casting is required.
    SetNthInput(0, const_cast<TInput *>( domain ));   
  }*/
 
 /***
 *	Generate Output Information
 **/
  template <class TInput, class TOutput>
  void Tractography<TInput,TOutput>
  ::GenerateOutputInformation()
  {

   // copy output information from input image
   Superclass::GenerateOutputInformation();

   // use user-specified output information
   //if ( this->GetInput() == NULL || m_OverrideOutputInformation ) {
	 OutputImagePointer		outputPtr = this->GetOutput();
     InputImageConstPointer	inputPtr = this->GetInput();
     
	 outputPtr->SetLargestPossibleRegion( inputPtr->GetLargestPossibleRegion() );
     outputPtr->SetVectorLength(m_FibersPerVoxel);
	 outputPtr->Allocate();
  }

	/*
		Compute fibers in the volume.
	*/
template <class TInput, class TOutput> 
void 
Tractography<TInput, TOutput>::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId){
		
	//std::cout<<"entro en tractografia"<<std::endl;
	InputImageConstPointer	inputPtr  = this->GetInput();
	OutputImagePointer		outputPtr = this->GetOutput();

	itk::VariableLengthVector<unsigned long> vzeros;
	vzeros.SetSize(m_FibersPerVoxel);
	vzeros.Fill(0);
	
	VariableLengthVector<unsigned long> stream_ant;
	stream_ant.SetSize(m_FibersPerVoxel);
		
	typedef itk::ImageRegionConstIteratorWithIndex< InputImageType>			InputIteratorType;
	InputIteratorType													it_in(inputPtr, outputRegionForThread);
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType>		OutputIteratorType;
	OutputIteratorType													it_out(outputPtr, outputRegionForThread);
	
	for (it_out.GoToBegin();!it_out.IsAtEnd();++it_out) {
		it_out.Set(vzeros);
	}
					
	StreamlineType streamline;
	float fa;
	InputIndexType index;
	
	m_StreamlineList.reserve(0);
	std::vector<StreamlineType> vectoraux;
	InputSpacingType spacing = inputPtr->GetSpacing();
	InputPointType	origin   = inputPtr->GetOrigin();
 	if (m_Connection == false) {
	  unsigned long stream_index=threadId*100000;
		std::cout << "Tractography brute force threadId " << threadId << std::endl;	
	  for (it_in.GoToBegin(),	it_out.GoToBegin();!it_in.IsAtEnd();++it_in,++it_out) {
		if((it_out.Get())[0]==0){
			fa=it_in.Get().GetFractionalAnisotropy();
			if(fa>=m_FaThreshold){
				streamline.clear();
				RungeKuttaIntegration(it_in.GetIndex(), streamline );				
				
				if(streamline.size()>=m_MinimumFiberLength){
					stream_index++;
				//	std::cout<<"fibra: "<<stream_index<<std::endl;
					vectoraux.push_back(streamline);
					for(unsigned int i=0; i<streamline.size(); ++i){
						index[0]=round((streamline[i][0]-origin[0])/spacing[0]);
						index[1]=round((streamline[i][1]-origin[1])/spacing[1]);
						index[2]=round((streamline[i][2]-origin[2])/spacing[2]);
						
						
						stream_ant=outputPtr->GetPixel(index);
						for(unsigned int j=0; j<stream_ant.Size(); ++j){
							if(stream_ant[j]==0){
								stream_ant[j]=stream_index;
								break;
							}
							if(stream_ant[j]==stream_index){
								break;
							}
						}
					
						outputPtr->SetPixel(index, stream_ant);
					}				
				}
			}
		}
	/*	if (stream_index>1000){
			break;
		}*/
	  }
	} else {
	  typedef itk::ImageRegionConstIterator< CharImageType>		CharIteratorType;
	  CharIteratorType													it_seeds(m_SeedsImage, inputPtr->GetLargestPossibleRegion());
	  unsigned long stream_index=threadId*100000;
		std::cout << "Tractography with conection threadId" << threadId << std::endl;	
	  for (it_in.GoToBegin(),	it_out.GoToBegin(), it_seeds.GoToBegin();!it_in.IsAtEnd();++it_in,++it_out, ++it_seeds) {
		if( it_seeds.Get()!=0 ){
			fa=it_in.Get().GetFractionalAnisotropy();
			if(fa>=m_FaThreshold){
				streamline.clear();
				RungeKuttaIntegration(it_in.GetIndex(), streamline );				
				
				if(streamline.size()>=m_MinimumFiberLength){
					stream_index++;
					//	std::cout<<"fibra: "<<stream_index<<std::endl;
					vectoraux.push_back(streamline);
					for(unsigned int i=0; i<streamline.size(); ++i){
						index[0]=round((streamline[i][0]-origin[0])/spacing[0]);
						index[1]=round((streamline[i][1]-origin[1])/spacing[1]);
						index[2]=round((streamline[i][2]-origin[2])/spacing[2]);
						
						
						stream_ant=outputPtr->GetPixel(index);
						for(unsigned int j=0; j<stream_ant.Size(); ++j){
							if(stream_ant[j]==0){
								stream_ant[j]=stream_index;
								break;
							}
							if(stream_ant[j]==stream_index){
								break;
							}
						}
						
						outputPtr->SetPixel(index, stream_ant);
					}				
				}
			}
		}
	  }
	}
	std::cout << "Done" << std::endl;	
	m_StreamlineList=vectoraux;
	//std::cout << "m_StreamlineList.size: " << m_StreamlineList.size() << std::endl;
};
	
template <class TInput, class TOutput> 
void 
Tractography<TInput, TOutput>::RungeKuttaIntegrationWithConnection( ){
		
		//std::cout<< "entro en RungeKuttaIntegrationWithConnection"<<std::endl;
		InputImageConstPointer	inputPtr  = this->GetInput();
	    InputSpacingType spacing = inputPtr->GetSpacing();
	    InputPointType	origin   = inputPtr->GetOrigin();
		OutputImagePointer		outputPtr = this->GetOutput();
	
	    outputPtr->SetLargestPossibleRegion( inputPtr->GetLargestPossibleRegion() );
	    outputPtr->SetVectorLength(m_FibersPerVoxel);
	    outputPtr->Allocate();
	
		itk::VariableLengthVector<unsigned long> vzeros;
		vzeros.SetSize(m_FibersPerVoxel);
		vzeros.Fill(0);
		
		VariableLengthVector<unsigned long> stream_ant;
		stream_ant.SetSize(m_FibersPerVoxel);

		typedef itk::ImageRegionConstIteratorWithIndex< InputImageType>			InputIteratorType;
		InputIteratorType													it_in(inputPtr, inputPtr->GetLargestPossibleRegion());
		typedef itk::ImageRegionIteratorWithIndex< OutputImageType>		OutputIteratorType;
		OutputIteratorType													it_out(outputPtr, inputPtr->GetLargestPossibleRegion());
    	typedef itk::ImageRegionConstIterator< CharImageType>		CharIteratorType;
	    CharIteratorType													it_seeds(m_SeedsImage, inputPtr->GetLargestPossibleRegion());

		for (it_out.GoToBegin();!it_out.IsAtEnd();++it_out) {
			it_out.Set(vzeros);
		}
		
		StreamlineType streamline;
		float fa;
		InputIndexType index;
		
		m_StreamlineList.reserve(0);
		std::vector<StreamlineType> vectoraux;
		
		
		unsigned long stream_index=0;
	    int count=0;
		for (it_in.GoToBegin(),	it_out.GoToBegin(), it_seeds.GoToBegin();!it_in.IsAtEnd();++it_in,++it_out, ++it_seeds) {
			if( it_seeds.Get()!=0 ){
				fa=it_in.Get().GetFractionalAnisotropy();
				if(fa>=m_FaThreshold){
					streamline.clear();
					RungeKuttaIntegration(it_in.GetIndex(), streamline );				
					
					if(streamline.size()>=m_MinimumFiberLength){
						stream_index++;
						//	std::cout<<"fibra: "<<stream_index<<std::endl;
						vectoraux.push_back(streamline);
						for(unsigned int i=0; i<streamline.size(); ++i){
							index[0]=round((streamline[i][0]-origin[0])/spacing[0]);
							index[1]=round((streamline[i][1]-origin[1])/spacing[1]);
							index[2]=round((streamline[i][2]-origin[2])/spacing[2]);

							stream_ant=outputPtr->GetPixel(index);
							for(unsigned int j=0; j<m_FibersPerVoxel; ++j){
								if(stream_ant[j]==0){
									stream_ant[j]=stream_index;
									break;
								}
								if(stream_ant[j]==stream_index){
									break;
								}
							}
							
							outputPtr->SetPixel(index, stream_ant);
							//std::cout << "setting in index:" << index << "vector.Size " << stream_ant.Size() << " vector[0]: "<< stream_ant[0] << std::endl;
						}				
					}
				}
			}
			/*	if (stream_index>1000){
			 break;
			 }*/
		}
		m_StreamlineList=vectoraux;
		std::cout << "m_StreamlineList.size: " << m_StreamlineList.size() << std::endl;
};
	
/**
*		Runge-Kutta Integration Method from a seed
**/
template <class TInput, class TOutput> 
void
Tractography<TInput, TOutput>::RungeKuttaIntegration(InputIndexType seed, StreamlineType& streamline){

	InputImageConstPointer inputPtr= this->GetInput();
	float eps=1e-06;

	InputSpacingType spacing = inputPtr->GetSpacing();
	InputPointType	origin   = inputPtr->GetOrigin();

	typename TInput::SizeType Size = inputPtr->GetLargestPossibleRegion().GetSize();
	
	InputPointType point;
	StreamlineType streamlineForward;
	/*StreamlineType streamlineBack;*/
//	StreamlineType streamline;

	float ki[4][3];
	float dist;
	float ksign;
	unsigned long voxelNumber=0;

	float curvature=0;
	float fa, fa_aux;
	float mod_vn;
	float mod_vector_ant;				
				
	InputIndexType voxelIndex;
	InputPointType rn;
  
	typedef itk::Vector<float,3>                       VectorType;
	
	typedef itk::DTITensor<float>::EigenValuesArrayType		EigenValuesType;
	EigenValuesType			eigval;
	itk::DTITensor<float>::RealType				eigvec[3];
	
	VectorType ern, vn, e_kp, e_k, e_neigh, rh, h, vector_ant;
	//float m_StepLength=0.1;
	h[0]=m_StepLength;
	h[1]=m_StepLength;
	h[2]=m_StepLength;
	
	if(m_CurvatureThreshold > 3.14159){
		m_CurvatureThreshold=(m_CurvatureThreshold/180)*3.14159;
	}
	
	typedef std::vector<InputIndexType>             IndexVectorType;
	IndexVectorType neigh;
	InputIndexType aux_neigh;

	fa=inputPtr->GetPixel(seed).GetFractionalAnisotropy();

	if(fa >= m_FaThreshold){
		voxelIndex=seed;
		if( curvature <= m_CurvatureThreshold){
			inputPtr->GetPixel(seed).ComputeEigenValues(eigval );
			inputPtr->GetPixel(seed).ComputeEigenVector(eigval[0], eigvec);
			ern[0]=eigvec[0];
			ern[1]=eigvec[1];
			ern[2]=eigvec[2];
			curvature=0;
		}
		for(unsigned s=0; s<2; s++){
			voxelIndex=seed;
			if(s==1){
				ern[0]=-ern[0];
				ern[1]=-ern[1];
				ern[2]=-ern[2];
				
				while(!streamlineForward.empty()){	
					streamline.push_back(streamlineForward.back());
					streamlineForward.pop_back();
				}
				streamline.pop_back();
			}
			vn[0]=ern[0];
			vn[1]=ern[1];
			vn[2]=ern[2];
			
			rn[0]=(seed[0])*spacing[0]+origin[0];
			rn[1]=(seed[1])*spacing[1]+origin[1];
			rn[2]=(seed[2])*spacing[2]+origin[2];
			
			fa=inputPtr->GetPixel(seed).GetFractionalAnisotropy();
			curvature=0;
	
			while(fa>=m_FaThreshold && curvature <= m_CurvatureThreshold && voxelNumber<100000){
				point=rn;
				voxelNumber++;
				if(s==0){
					streamlineForward.push_back(point);
				}else{
					streamline.push_back(point);
				}
				
				//Buscar una forma mejor de buscar los pixeles que lo rodean.
				aux_neigh[0]=floor((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
				
				aux_neigh[0]=ceil((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				e_k.Fill(0);	
				for(unsigned n=0; n<8; n++){
					if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
						dist=sqrt((rn[0]-(spacing[0]*neigh[n][0]+origin[0]))*(rn[0]-(spacing[0]*neigh[n][0]+origin[0]))+(rn[1]-(spacing[1]*neigh[n][1]+origin[1]))*(rn[1]-(spacing[1]*neigh[n][1]+origin[1]))+(rn[2]-(spacing[2]*neigh[n][2]+origin[2]))*(rn[2]-(spacing[2]*neigh[n][2]+origin[2])));
						fa_aux=inputPtr->GetPixel(neigh[n]).GetFractionalAnisotropy();
						if(fa_aux>m_FaThreshold){
							inputPtr->GetPixel(neigh[n]).ComputeEigenValues(eigval );
							inputPtr->GetPixel(neigh[n]).ComputeEigenVector(eigval[0], eigvec);
							e_neigh[0]=eigvec[0];
							e_neigh[1]=eigvec[1];
							e_neigh[2]=eigvec[2];
							if(dist<eps){
								e_k[0]=eigvec[0];
								e_k[1]=eigvec[1];
								e_k[2]=eigvec[2];
								break;					
							}
							e_k=e_k+fa_aux*(1/dist)*e_neigh;
						}
					}
				}
				e_k.Normalize();
				neigh.clear();
				
				inputPtr->GetPixel(voxelIndex).ComputeEigenValues(eigval );
				inputPtr->GetPixel(voxelIndex).ComputeEigenVector(eigval[0], eigvec);
				
				e_kp[0]=eigvec[0];
				e_kp[1]=eigvec[1];
				e_kp[2]=eigvec[2];
				e_k=e_kp;
			
				//EMMA->Cambio la forma de calcular los valores de k;
				//	for(unsigned p=0; p<3; p++){
				//ksign=vn*ern;
				ksign=vn*e_k;
				if(fabs(ksign)>eps){
					ksign=(ksign/fabs(ksign));
				}else{
					ksign=1;
				}
				
				for(unsigned n=0;n<3;n++){
					ki[0][n]=ksign*e_k[n];
					rh[n]=point[n]+((h[n]/2)*ki[0][n]);
				}
								
				//Buscar una forma mejor de buscar los pixeles que lo rodean.
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
				
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
			
				e_k[0]=0;
				e_k[1]=0;
				e_k[2]=0;
				
				for(unsigned n=0; n<8; n++){
					if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
						dist=sqrt((rh[0]-(spacing[0]*neigh[n][0]+origin[0]))*(rh[0]-(spacing[0]*neigh[n][0]+origin[0]))+(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))*(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))+(rh[2]-(spacing[2]*neigh[n][2]+origin[2]))*(rh[2]-(spacing[2]*neigh[n][2]+origin[2])));
						fa_aux=inputPtr->GetPixel(neigh[n]).GetFractionalAnisotropy();
						if(fa_aux>m_FaThreshold){
							inputPtr->GetPixel(neigh[n]).ComputeEigenValues(eigval );
							inputPtr->GetPixel(neigh[n]).ComputeEigenVector(eigval[0], eigvec);
							e_neigh[0]=eigvec[0];
							e_neigh[1]=eigvec[1];
							e_neigh[2]=eigvec[2];
							if(dist<eps){
								e_k[0]=eigvec[0];
								e_k[1]=eigvec[1];
								e_k[2]=eigvec[2];
								break;					
							}
							e_k=e_k+fa_aux*(1/dist)*e_neigh;
						}
					}
				}
				e_k.Normalize();
				neigh.clear();
				
				ksign=vn*e_k;
				if(fabs(ksign)>eps){
					ksign=(ksign/fabs(ksign));
				}else{
					ksign=1;
				}
				
				for(unsigned n=0;n<3;n++){
					ki[1][n]=ksign*e_k[n];
					rh[n]=point[n]+((h[n]/2)*ki[1][n]);
				}
	
				//Buscar una forma mejor de buscar los pixeles que lo rodean.
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
				
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
						
				e_k.Fill(0);		
				for(unsigned n=0; n<8; n++){
					if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
						dist=sqrt((rh[0]-(spacing[0]*neigh[n][0]+origin[0]))*(rh[0]-(spacing[0]*neigh[n][0]+origin[0]))+(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))*(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))+(rh[2]-(spacing[2]*neigh[n][2]+origin[2]))*(rh[2]-(spacing[2]*neigh[n][2]+origin[2])));
						fa_aux=inputPtr->GetPixel(neigh[n]).GetFractionalAnisotropy();
						if(fa_aux>m_FaThreshold){

							inputPtr->GetPixel(neigh[n]).ComputeEigenValues(eigval );
							inputPtr->GetPixel(neigh[n]).ComputeEigenVector(eigval[0], eigvec);
					
							e_neigh[0]=eigvec[0];	
							e_neigh[1]=eigvec[1];
							e_neigh[2]=eigvec[2];
							if(dist<eps){
								e_k[0]=eigvec[0];
								e_k[1]=eigvec[1];
								e_k[2]=eigvec[2];
								break;					
							}
							e_k=e_k+fa_aux*(1/dist)*e_neigh;
						}
					}
				}
				e_k.Normalize();
				neigh.clear();
				
			
				ksign=vn*e_k;
				if(fabs(ksign)>eps){
					ksign=(ksign/fabs(ksign));
				}else{
					ksign=1;
				}
				
				for(unsigned n=0;n<3;n++){
					ki[2][n]=ksign*e_k[n];
					rh[n]=point[n]+(h[n]*ki[2][n]);
				}
	
				
				//Buscar una forma mejor de buscar los pixeles que lo rodean.
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
				
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				e_k.Fill(0);	
				for(unsigned n=0; n<8; n++){
					if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
						dist=sqrt((rh[0]-(spacing[0]*neigh[n][0]+origin[0]))*(rh[0]-(spacing[0]*neigh[n][0]+origin[0]))+(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))*(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))+(rh[2]-(spacing[2]*neigh[n][2]+origin[2]))*(rh[2]-(spacing[2]*neigh[n][2]+origin[2])));
					
						fa_aux=inputPtr->GetPixel(neigh[n]).GetFractionalAnisotropy();
						if(fa_aux>m_FaThreshold){
							inputPtr->GetPixel(neigh[n]).ComputeEigenValues(eigval );
							inputPtr->GetPixel(neigh[n]).ComputeEigenVector(eigval[0], eigvec);
							e_neigh[0]=eigvec[0];
							e_neigh[1]=eigvec[1];
							e_neigh[2]=eigvec[2];
							if(dist<eps){
								e_k[0]=eigvec[0];
								e_k[1]=eigvec[1];
								e_k[2]=eigvec[2];		
								break;			
							}
							e_k=e_k+fa_aux*(1/dist)*e_neigh;
						}
					}
				}
				e_k.Normalize();
				neigh.clear();
				
				ksign=vn*e_k;
				if(fabs(ksign)>eps){
					ksign=(ksign/fabs(ksign));
				}else{
					ksign=1;
				}
				
				for(unsigned n=0;n<3;n++){
					ki[3][n]=ksign*e_k[n];
				}	
				// Cambio hecho por Ruben el 3-6-08
				//vector_ant=ern;
				vector_ant=vn;
			
				for(unsigned n=0; n<3; n++){
					vn[n]=(ki[0][n]+2*ki[1][n]+2*ki[2][n]+ki[3][n])/6;
					rn[n]=rn[n]+h[n]*vn[n];
				}
			
				if((vn*vector_ant)<0){
					vn[0]=-vn[0];
					vn[1]=-vn[1];
					vn[2]=-vn[2];
				}
			
				for(unsigned n=0; n<3; n++){
				//	vn[n]=(ki[0][n]+2*ki[1][n]+2*ki[2][n]+ki[3][n])/6;
					rn[n]=rn[n]+h[n]*vn[n];
				}
				
				
				
				
				voxelIndex[0]=round((rn[0]-origin[0])/spacing[0]);
				voxelIndex[1]=round((rn[1]-origin[1])/spacing[1]);
				voxelIndex[2]=round((rn[2]-origin[2])/spacing[2]);
			    if (voxelIndex[0] < Size[0] && 
				    voxelIndex[1] < Size[1] && 
					voxelIndex[2] < Size[2] &&
					voxelIndex[0] >  0  && 
				    voxelIndex[1] >  0 && 
					voxelIndex[2] >	 0 ) {
				    fa=inputPtr->GetPixel(voxelIndex).GetFractionalAnisotropy(); 
				} else {
				  fa = 0;
				} 
				
				curvature=fabs(vn*vector_ant);
				mod_vn=sqrt(vn*vn);
				mod_vector_ant=sqrt(vector_ant*vector_ant);				
				curvature=curvature/(mod_vn*mod_vector_ant);
				
				if(curvature>1){
					curvature=1;
				};
				curvature=fabs(acos( curvature));
				
			}
		}
	}
};
	
template <class TInput, class TOutput> 
void
Tractography<TInput, TOutput>::RungeKuttaIntegrationNew(InputIndexType seed, StreamlineType& streamline){

	InputImageConstPointer inputPtr= this->GetInput();
	float eps=1e-06;

	InputSpacingType spacing = inputPtr->GetSpacing();
	InputPointType	origin   = inputPtr->GetOrigin();

	typename TInput::SizeType Size = inputPtr->GetLargestPossibleRegion().GetSize();
	
	InputPointType point;
	StreamlineType streamlineForward;
	/*StreamlineType streamlineBack;*/
//	StreamlineType streamline;

	float ki[4][3];
	float dist;
	float ksign;
	unsigned long voxelNumber=0;

	float curvature=0;
	float fa, fa_aux;
	float mod_vn;
	float mod_vector_ant;				
				
	InputIndexType voxelIndex;
	InputPointType rn;
  
	typedef itk::Vector<float,3>                       VectorType;
	
	typedef itk::DTITensor<float>::EigenValuesArrayType		EigenValuesType;
	EigenValuesType			eigval;
	//itk::DTITensor<float>::RealType				eigvec[3];
	itk::SymmetricSecondRankTensor<float,3> symTensor;
	itk::SymmetricSecondRankTensor<float,3>::EigenValuesArrayType eigVal;
	itk::SymmetricSecondRankTensor<float,3>::EigenVectorsMatrixType eigVec;
	
	VectorType ern, vn, e_kp, e_k, e_neigh, rh, h, vector_ant;
	//float m_StepLength=0.1;
	h[0]=m_StepLength;
	h[1]=m_StepLength;
	h[2]=m_StepLength;
	
	if(m_CurvatureThreshold > 3.14159){
		m_CurvatureThreshold=(m_CurvatureThreshold/180)*3.14159;
	}
	
	typedef std::vector<InputIndexType>             IndexVectorType;
	IndexVectorType neigh;
	InputIndexType aux_neigh;

	fa=inputPtr->GetPixel(seed).GetFractionalAnisotropy();

	if(fa >= m_FaThreshold){
		voxelIndex=seed;
		if( curvature <= m_CurvatureThreshold){
			//inputPtr->GetPixel(seed).ComputeEigenValues(eigval );
			//inputPtr->GetPixel(seed).ComputeEigenVector(eigval[0], eigvec);
			;
			for(unsigned i=0; i<6; i++){
   		       symTensor[i]=inputPtr->GetPixel(seed)[i];
	        }
			symTensor.ComputeEigenAnalysis(eigVal, eigVec);
			ern[0]=eigVec[2][0];
			ern[1]=eigVec[2][1];
			ern[2]=eigVec[2][2];
			curvature=0;
		}
		for(unsigned s=0; s<2; s++){
			voxelIndex=seed;
			if(s==1){
				ern[0]=-ern[0];
				ern[1]=-ern[1];
				ern[2]=-ern[2];
				
				while(!streamlineForward.empty()){	
					streamline.push_back(streamlineForward.back());
					streamlineForward.pop_back();
				}
				streamline.pop_back();
			}
			vn[0]=ern[0];
			vn[1]=ern[1];
			vn[2]=ern[2];
			
			rn[0]=(seed[0])*spacing[0]+origin[0];
			rn[1]=(seed[1])*spacing[1]+origin[1];
			rn[2]=(seed[2])*spacing[2]+origin[2];
			
			fa=inputPtr->GetPixel(seed).GetFractionalAnisotropy();
			curvature=0;
	
			while(fa>=m_FaThreshold && curvature <= m_CurvatureThreshold && voxelNumber<100000){
				point=rn;
				voxelNumber++;
				if(s==0){
					streamlineForward.push_back(point);
				}else{
					streamline.push_back(point);
				}
				
				//Buscar una forma mejor de buscar los pixeles que lo rodean.
				aux_neigh[0]=floor((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
				
				aux_neigh[0]=ceil((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rn[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rn[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rn[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				e_k.Fill(0);	
				for(unsigned n=0; n<8; n++){
					if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
						dist=sqrt((rn[0]-(spacing[0]*neigh[n][0]+origin[0]))*(rn[0]-(spacing[0]*neigh[n][0]+origin[0]))+(rn[1]-(spacing[1]*neigh[n][1]+origin[1]))*(rn[1]-(spacing[1]*neigh[n][1]+origin[1]))+(rn[2]-(spacing[2]*neigh[n][2]+origin[2]))*(rn[2]-(spacing[2]*neigh[n][2]+origin[2])));
						fa_aux=inputPtr->GetPixel(neigh[n]).GetFractionalAnisotropy();
						if(fa_aux>m_FaThreshold){
						
							//inputPtr->GetPixel(neigh[n]).ComputeEigenValues(eigval );
							//inputPtr->GetPixel(neigh[n]).ComputeEigenVector(eigval[0], eigvec);
							for(unsigned i=0; i<6; i++){
   		                      symTensor[i]=inputPtr->GetPixel(neigh[n])[i];
	                        }
			                symTensor.ComputeEigenAnalysis(eigVal, eigVec);
							e_neigh[0]=eigVec[2][0];
							e_neigh[1]=eigVec[2][1];
							e_neigh[2]=eigVec[2][2];
							if(dist<eps){
								e_k[0]=eigVec[2][0];
								e_k[1]=eigVec[2][1];
								e_k[2]=eigVec[2][2];
								break;					
							}
							e_k=e_k+fa_aux*(1/dist)*e_neigh;
						}
					}
				}
				e_k.Normalize();
				neigh.clear();
				
				//inputPtr->GetPixel(voxelIndex).ComputeEigenValues(eigval );
				//inputPtr->GetPixel(voxelIndex).ComputeEigenVector(eigval[0], eigvec);
				for(unsigned i=0; i<6; i++){
					symTensor[i]=inputPtr->GetPixel(voxelIndex)[i];
				}
				symTensor.ComputeEigenAnalysis(eigVal, eigVec);
				
							
				e_kp[0]=eigVec[2][0];
				e_kp[1]=eigVec[2][1];
				e_kp[2]=eigVec[2][2];
				e_k=e_kp;
			
				//EMMA->Cambio la forma de calcular los valores de k;
				//	for(unsigned p=0; p<3; p++){
				//ksign=vn*ern;
				ksign=vn*e_k;
				if(fabs(ksign)>eps){
					ksign=(ksign/fabs(ksign));
				}else{
					ksign=1;
				}
				
				for(unsigned n=0;n<3;n++){
					ki[0][n]=ksign*e_k[n];
					rh[n]=point[n]+((h[n]/2)*ki[0][n]);
				}
								
				//Buscar una forma mejor de buscar los pixeles que lo rodean.
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
				
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
			
				e_k[0]=0;
				e_k[1]=0;
				e_k[2]=0;
				
				for(unsigned n=0; n<8; n++){
					if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
						dist=sqrt((rh[0]-(spacing[0]*neigh[n][0]+origin[0]))*(rh[0]-(spacing[0]*neigh[n][0]+origin[0]))+(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))*(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))+(rh[2]-(spacing[2]*neigh[n][2]+origin[2]))*(rh[2]-(spacing[2]*neigh[n][2]+origin[2])));
						fa_aux=inputPtr->GetPixel(neigh[n]).GetFractionalAnisotropy();
						if(fa_aux>m_FaThreshold){
							//inputPtr->GetPixel(neigh[n]).ComputeEigenValues(eigval );
							//inputPtr->GetPixel(neigh[n]).ComputeEigenVector(eigval[0], eigvec);
							for(unsigned i=0; i<6; i++){
   		                      symTensor[i]=inputPtr->GetPixel(neigh[n])[i];
	                        }
			                symTensor.ComputeEigenAnalysis(eigVal, eigVec);
							e_neigh[0]=eigVec[2][0];
							e_neigh[1]=eigVec[2][1];
							e_neigh[2]=eigVec[2][2];
							if(dist<eps){
								e_k[0]=eigVec[2][0];
								e_k[1]=eigVec[2][1];
								e_k[2]=eigVec[2][2];
								break;					
							}
							e_k=e_k+fa_aux*(1/dist)*e_neigh;
						}
					}
				}
				e_k.Normalize();
				neigh.clear();
				
				ksign=vn*e_k;
				if(fabs(ksign)>eps){
					ksign=(ksign/fabs(ksign));
				}else{
					ksign=1;
				}
				
				for(unsigned n=0;n<3;n++){
					ki[1][n]=ksign*e_k[n];
					rh[n]=point[n]+((h[n]/2)*ki[1][n]);
				}
	
				//Buscar una forma mejor de buscar los pixeles que lo rodean.
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
				
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
						
				e_k.Fill(0);		
				for(unsigned n=0; n<8; n++){
					if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
						dist=sqrt((rh[0]-(spacing[0]*neigh[n][0]+origin[0]))*(rh[0]-(spacing[0]*neigh[n][0]+origin[0]))+(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))*(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))+(rh[2]-(spacing[2]*neigh[n][2]+origin[2]))*(rh[2]-(spacing[2]*neigh[n][2]+origin[2])));
						fa_aux=inputPtr->GetPixel(neigh[n]).GetFractionalAnisotropy();
						if(fa_aux>m_FaThreshold){

							//inputPtr->GetPixel(neigh[n]).ComputeEigenValues(eigval );
							//inputPtr->GetPixel(neigh[n]).ComputeEigenVector(eigval[0], eigvec);
					        for(unsigned i=0; i<6; i++){
   		                      symTensor[i]=inputPtr->GetPixel(neigh[n])[i];
	                        }
			                symTensor.ComputeEigenAnalysis(eigVal, eigVec);
							e_neigh[0]=eigVec[2][0];	
							e_neigh[1]=eigVec[2][1];
							e_neigh[2]=eigVec[2][2];
							if(dist<eps){
								e_k[0]=eigVec[2][0];
								e_k[1]=eigVec[2][1];
								e_k[2]=eigVec[2][2];
								break;					
							}
							e_k=e_k+fa_aux*(1/dist)*e_neigh;
						}
					}
				}
				e_k.Normalize();
				neigh.clear();
				
			
				ksign=vn*e_k;
				if(fabs(ksign)>eps){
					ksign=(ksign/fabs(ksign));
				}else{
					ksign=1;
				}
				
				for(unsigned n=0;n<3;n++){
					ki[2][n]=ksign*e_k[n];
					rh[n]=point[n]+(h[n]*ki[2][n]);
				}
	
				
				//Buscar una forma mejor de buscar los pixeles que lo rodean.
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=floor((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
				
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=floor((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=floor((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				aux_neigh[0]=ceil((rh[0]-origin[0])/spacing[0]);
				aux_neigh[1]=ceil((rh[1]-origin[1])/spacing[1]);
				aux_neigh[2]=ceil((rh[2]-origin[2])/spacing[2]);
				neigh.push_back(aux_neigh);
					
				e_k.Fill(0);	
				for(unsigned n=0; n<8; n++){
					if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
						dist=sqrt((rh[0]-(spacing[0]*neigh[n][0]+origin[0]))*(rh[0]-(spacing[0]*neigh[n][0]+origin[0]))+(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))*(rh[1]-(spacing[1]*neigh[n][1]+origin[1]))+(rh[2]-(spacing[2]*neigh[n][2]+origin[2]))*(rh[2]-(spacing[2]*neigh[n][2]+origin[2])));
					
						fa_aux=inputPtr->GetPixel(neigh[n]).GetFractionalAnisotropy();
						if(fa_aux>m_FaThreshold){
							//inputPtr->GetPixel(neigh[n]).ComputeEigenValues(eigval );
							//inputPtr->GetPixel(neigh[n]).ComputeEigenVector(eigval[0], eigvec);
							for(unsigned i=0; i<6; i++){
   		                      symTensor[i]=inputPtr->GetPixel(neigh[n])[i];
	                        }
			                symTensor.ComputeEigenAnalysis(eigVal, eigVec);
							e_neigh[0]=eigVec[2][0];
							e_neigh[1]=eigVec[2][1];
							e_neigh[2]=eigVec[2][2];
							if(dist<eps){
								e_k[0]=eigVec[2][0];
								e_k[1]=eigVec[2][1];
								e_k[2]=eigVec[2][2];		
								break;			
							}
							e_k=e_k+fa_aux*(1/dist)*e_neigh;
						}
					}
				}
				e_k.Normalize();
				neigh.clear();
				
				ksign=vn*e_k;
				if(fabs(ksign)>eps){
					ksign=(ksign/fabs(ksign));
				}else{
					ksign=1;
				}
				
				for(unsigned n=0;n<3;n++){
					ki[3][n]=ksign*e_k[n];
				}	
				// Cambio hecho por Ruben el 3-6-08
				//vector_ant=ern;
				vector_ant=vn;
			
				for(unsigned n=0; n<3; n++){
					vn[n]=(ki[0][n]+2*ki[1][n]+2*ki[2][n]+ki[3][n])/6;
					rn[n]=rn[n]+h[n]*vn[n];
				}
			
				if((vn*vector_ant)<0){
					vn[0]=-vn[0];
					vn[1]=-vn[1];
					vn[2]=-vn[2];
				}
			
				for(unsigned n=0; n<3; n++){
				//	vn[n]=(ki[0][n]+2*ki[1][n]+2*ki[2][n]+ki[3][n])/6;
					rn[n]=rn[n]+h[n]*vn[n];
				}
				
				
				
				
				voxelIndex[0]=round((rn[0]-origin[0])/spacing[0]);
				voxelIndex[1]=round((rn[1]-origin[1])/spacing[1]);
				voxelIndex[2]=round((rn[2]-origin[2])/spacing[2]);
			    if (voxelIndex[0] < Size[0] && 
				    voxelIndex[1] < Size[1] && 
					voxelIndex[2] < Size[2] &&
					voxelIndex[0] >  0  && 
				    voxelIndex[1] >  0 && 
					voxelIndex[2] >	 0 ) {
				    fa=inputPtr->GetPixel(voxelIndex).GetFractionalAnisotropy(); 
				} else {
				  fa = 0;
				} 
				
				curvature=fabs(vn*vector_ant);
				mod_vn=sqrt(vn*vn);
				mod_vector_ant=sqrt(vector_ant*vector_ant);				
				curvature=curvature/(mod_vn*mod_vector_ant);
				
				if(curvature>1){
					curvature=1;
				};
				curvature=fabs(acos( curvature));
				
			}
		}
	}
};
	
	
 /* template <class TInput, class TOutput  >
  void Tractography<TInput,TOutput >
  ::GenerateData( ) {
  
  InputPointer m_Input = this->GetInput();
  OutputPointer m_Output = this->GetOutput(0);
  m_Output->SetBufferedRegion( m_Input->GetRequestedRegion() );
  m_Output->Allocate();
  
  RungeKuttaIntegration();

  
}*/

}
   

#endif

