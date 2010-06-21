/* Copyright (c) Ruben Cardenes Almeida 27/08/2002 */

#ifndef _geodesicPath3D_txx
#define _geodesicPath3D_txx

#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkVector.h>
#include <itkNumericTraits.h>
#include <itkImageFileWriter.h>
#include "geodesicPath3D.h"

#define PI 3.1415927
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

namespace itk {

   /**
   * Constructor
   */
  
  template <class TInput, class TOutput> 
  geodesicPath3D<TInput, TOutput> 
  ::geodesicPath3D()
    : m_Heap()
  {
   
    this->SetNumberOfRequiredInputs( 1 );
    this->SetNumberOfRequiredOutputs( 1 );
        
    /**
     * A pointer to the output image is created and conected to the output
     **/    
    OutputPointer output = OutputType::New();
    this->SetNthOutput( 0, output.GetPointer() );
	OutputSizeType outputSize;
    outputSize.Fill( 16 );
	m_OverrideOutputInformation = false;
	m_OutputRegion.SetSize( outputSize );
   
	pi = 3.1415927;
	m_max_num_points = 30;
	m_Verbose = false;
   }// End Constructor

   /**
   * SETINPUT
   **/
  template <class TInput, class TOutput  >
  void geodesicPath3D<TInput,TOutput >
  ::SetInput(const  TInput * domain ) 
  {
    // Process object is not const-correct so the const casting is required.
    SetNthInput(0, const_cast<TInput *>( domain ));   
  }

  template <class TInput, class TOutput>
  void geodesicPath3D<TInput,TOutput>
  ::GenerateOutputInformation()
  {

   // copy output information from input image
   Superclass::GenerateOutputInformation();

   // use user-specified output information
   if ( this->GetInput() == NULL || m_OverrideOutputInformation ) {
     OutputPointer output = this->GetOutput();
     output->SetLargestPossibleRegion( m_OutputRegion );
     output->SetSpacing( m_OutputSpacing );
     output->SetOrigin( m_OutputOrigin );
   }    
  }
  
  template <class TInput, class TOutput  >
  void geodesicPath3D<TInput,TOutput>
  ::SetSeeds(NodeContainer * points) 
  {
    m_seeds = points;
	this->Modified();
  }

 // template <class TInput, class TOutput  >
//  NodeContainer* geodesicPath3D<TInput,TOutput>
//  ::GetSeeds() 
//  {
//    return m_seeds_out;
//  }

  template <class TInput, class TOutput  >
  void geodesicPath3D<TInput,TOutput >
  ::GetNewCoordinates3D( int *newx,int *newy,int* newz,int x,int y,int z, 
                     float gradientx,float gradienty,float gradientz,
					 float angle_fi_error,float angle_theta_error,
					 float *new_fi_error,float *new_theta_error) {

  double theta,newtheta,fi,newfi,dy,dx;
  (*new_fi_error) = 0;
  (*new_theta_error) = 0;
  theta = acos((double)gradientz); // (debe salir un theta entre 0 y PI)
  newtheta = theta + angle_theta_error;
  
  if (newtheta <= PI/8) {
    (*newz) = z+1; 
	(*newx) = x;
	(*newy) = y;
	(*new_theta_error) = newtheta;
  }
  if (newtheta <= 3*PI/8 && newtheta > PI/8) {
    (*newz) = z+1; 
	(*new_theta_error) = newtheta - PI/4;
  }
  if (newtheta > 3*PI/8 && newtheta < 5*PI/8) {
    (*newz) = z; 
    (*new_theta_error) = newtheta - PI/2;
  }
  if (newtheta >=5*PI/8 && newtheta < 7*PI/8) {
    (*newz) = z-1; 
    (*new_theta_error) = newtheta - 3*PI/4;
  }
  if (newtheta >=7*PI/8) {
    (*newz) = z-1; 
	(*newx) = x;
	(*newy) = y;
    (*new_theta_error) = newtheta - PI;
  }
  
  if (fabs(sin(theta)) < (double)fabs(gradientx)) {
    if (gradientx > 0) {
      fi = 0;
    }
    if (gradientx < 0) {
      fi = PI;
	}
  } else { 
    if (theta != 0 && theta != PI) {
      fi = acos((double)gradientx/sin(theta));
	}
  }
  if (gradienty < 0) {
	fi = 2 * PI - fi;
  }
  
  //if (gradientx != 0) {
//	fi = atan((double)gradienty/gradientx);
//  } else {
//	fi = PI/2;
//  }
//   
//  if (gradienty < 0) {
//	fi = 2 * PI - fi;
//  }
  
  newfi = fi + angle_fi_error;
  dy = sin(newfi);
  dx = cos(newfi);
  
  if (newtheta < PI/8 || newtheta > 7*PI/8) {
     // Todo decidido. Calcular error en fi. 
	 (*new_fi_error) = 0;
  } else {
    if (dy >=0) { /* theta -> (0,PI)*/
      if (fabs(dx) <= 0.3826) { /*theta -> (3PI/8,5PI/8) */
        (*newx) = x; 
        (*newy) = y-1; 
        (*new_fi_error) = newfi - PI / 2;
      } 
      if (dx >= 0.9238) { /*theta -> (0,PI/8) */
        (*newx) = x+1;
        (*newy) = y;
        (*new_fi_error) = newfi;
      }
      if (dx <= -0.9238) { /*theta -> (7PI/8,PI) */
        (*newx) = x-1;
        (*newy) = y;
        (*new_fi_error) = newfi - PI;
      }
      if (dx <= 0.9238 && dx >= 0.3826) { /*theta -> (PI/8,3PI/8) */
        (*newx) = x+1;
        (*newy) = y-1;
        (*new_fi_error) = newfi - PI / 4;
      }
      if (dx >= -0.9238 && dx <= -0.3826) { /*theta -> (5PI/8,7PI/8) */
        (*newx) = x-1;
        (*newy) = y-1;
        (*new_fi_error) = newfi - 3*PI / 4;
      }
    }
    if (dy <0) {
      if (fabs(dx) <= 0.3826) {
        (*newx) = x;
        (*newy) = y+1;
        (*new_fi_error) = newfi - 3*PI / 2;
      }
      if (dx >= 0.9238) {
        (*newx) = x+1;
        (*newy) = y;
        (*new_fi_error) = newfi;
      }
      if (dx <= -0.9238) {
        (*newx) = x-1;
        (*newy) = y;
        (*new_fi_error) = newfi - PI;
      }
      if (dx <= 0.9238 && dx >= 0.3826) {
        (*newx) = x+1;
        (*newy) = y+1;
        (*new_fi_error) = newfi -  7*PI / 4;
      }
      if (dx >= -0.9238 && dx <= -0.3826) {
        (*newx) = x-1;
        (*newy) = y+1; 
        (*new_fi_error) = newfi - 5*PI / 4;
      }
    }
  }
  
  return;
  }
  
  //int geodesicPath3D(int max1, int max2,int max3,float *maps,struct element Proto_end, char* proto, float ***gradientx, float*** gradienty,float*** gradientz) {
  template <class TInput, class TOutput  >
  void geodesicPath3D<TInput,TOutput >
  ::GenerateData( ) {
  
  InputPointer input = this->GetInput();
  OutputPointer output = this->GetOutput(0);
  output->SetBufferedRegion( input->GetRequestedRegion() );
  output->Allocate();

  InputPixelType dist_act;
  InputIndexType mapindex,mapindex_min,newmapindex;
  int count;

  NodeType node,node_new;

  InputSpacingType spacing = input->GetSpacing();
  InputPointType	origin   = input->GetOrigin();
  InputPointType rn;
  
  typedef typename itk::NeighborhoodIterator<InputType> inputNeighborIteratorType;
  typedef typename inputNeighborIteratorType::RadiusType radiusType;
  radiusType radius;
  radius.Fill(1);
  inputNeighborIteratorType inputNIt = inputNeighborIteratorType(radius, output, input->GetRequestedRegion());
  //typedef typename outputConstNeighborIteratorType::OffsetType offsetType;
  //offsetType offset;
  //offsetType v_offset[Dimension*2];
  
  m_PointsOut.clear();
  m_PointsOutFloat.clear();
  
  std::vector<OutputIndexType> PointsOutElement;
  PointsOutElement.clear();
  std::vector<InputPointType> PointsOutElementFloat;
  PointsOutElementFloat.clear();
  if (m_Metodo == 1) {
    if ( m_seeds )  {
	  typename NodeContainer::ConstIterator pointsIter = m_seeds->Begin();
	  typename NodeContainer::ConstIterator pointsEnd = m_seeds->End();	
	  for ( ; pointsIter != pointsEnd; ++pointsIter ) {
	    PointsOutElement.clear();
	    node = pointsIter.Value();
	    mapindex = node.GetIndex();
        dist_act = input->GetPixel(mapindex);
	    //std::cout << "Incluida semilla en " << node.GetIndex() << " dist_act: " << dist_act << std::endl;
	    count = 0;
        while(dist_act > 0) {
	      inputNIt.SetLocation(mapindex);
	      //for (unsigned int k = 0; k<inputNIt.GetNeighborhood().Size(); k++) {
		  for (unsigned int k = 0; k<27; k++) {
		    newmapindex = inputNIt.GetIndex(k);
		    if (newmapindex != mapindex &&
		        input->GetPixel(newmapindex) < dist_act && 
		        input->GetPixel(newmapindex) > -1) {	
			  mapindex_min = newmapindex; 
			  dist_act = input->GetPixel(newmapindex);
			  std::cout << "Newmapindex: " << newmapindex << " Dist Act: " << dist_act << std::endl;
		    }	
	      }
          std::cout << "Dist Act min: " << dist_act << std::endl;
          if (mapindex == mapindex_min) {
            std::cout << "Error" << std::endl;
	        return;
          }
       
          mapindex = mapindex_min;
		  node_new.SetIndex(mapindex);
		  node_new.SetValue(node.GetValue());
		  PointsOutElement.push_back(mapindex);
          output->SetPixel(mapindex,node.GetValue());
          count++;
       } // end while dist_act
	   //std::cout << "Num elem in fiber " << m_PointsOut.size() << std::endl;
	   m_PointsOut.push_back(PointsOutElement);
	   //std::cout << "Num elem in path = " << count << std::endl;
	  } // end for pointsIter
    } // end if m_seeds 
  } // end if m_Metodo == 1

  if (m_Metodo == 2) {
    // Hacemos los gradientes:
    typename GradFilterType::Pointer gx = GradFilterType::New(); 
    gx->SetInput(input);
	gx->SetOrder(1);
	gx->SetDirection(0);
    gx->Update();
		
    typename GradFilterType::Pointer gy = GradFilterType::New(); 
    gy->SetInput(input);
	gy->SetOrder(1);
	gy->SetDirection(1);
    gy->Update();
		
	typename GradFilterType::Pointer gz = GradFilterType::New(); 
    gz->SetInput(input);
	gz->SetOrder(1);
	gz->SetDirection(2);
    gz->Update();
	
    VectorFloatType float_index;
    
	if ( m_seeds )  {
	  int newx,newy,newz;
      float new_fi_error=0,new_theta_error=0,angle_fi_error=0,angle_theta_error=0;
      float Gx,Gy,Gz,N;
	  typename NodeContainer::ConstIterator pointsIter = m_seeds->Begin();
	  typename NodeContainer::ConstIterator pointsEnd = m_seeds->End();
	  for ( ; pointsIter != pointsEnd; ++pointsIter ) {
	    PointsOutElement.clear();
		PointsOutElementFloat.clear();
	    node = pointsIter.Value();
	    std::cout << "Metodo 2. Incluida semilla en " << node.GetIndex() << std::endl;
	    mapindex = node.GetIndex();
        dist_act = input->GetPixel(mapindex);
	    count = 0;
		output->SetPixel(mapindex,2);
        while (dist_act > 0) {
		  Gx = gx->GetOutput()->GetPixel(mapindex);
		  Gy = gy->GetOutput()->GetPixel(mapindex);
		  Gz = gz->GetOutput()->GetPixel(mapindex);
		  //printf("Gx %f Gy %f Gz %f\n",Gx,Gy,Gz);
		  N = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
		  Gx = Gx/N;
		  Gy = Gy/N;
		  Gz = Gz/N;
		  float_index[0] = (float)mapindex[0]-Gx;
		  float_index[1] = (float)mapindex[1]+Gy;
		  float_index[2] = (float)mapindex[2]-Gz;
          GetNewCoordinates3D(&newx,&newy,&newz,mapindex[0],mapindex[1],mapindex[2],
	                     	  -Gx,Gy,-Gz,
							  angle_fi_error,angle_theta_error,&new_fi_error,&new_theta_error);
          newmapindex[0] = newx;
		  newmapindex[1] = newy;
		  newmapindex[2] = newz;
		  
	      //printf("newmapindex %d, newx %d, newy %d, newz %d\n", newmapindex, newx, newy, newz);
	      // Lo siguiente es porque no podemos salirnos del campo de distancias, si el nuevo pixel esta fuera 
	      // buscamos el vecino mas adecuado
          if ( input->GetPixel(newmapindex) < 0 ) {
		    // printf("Entrando en chungo!\n"); 
		    inputNIt.SetLocation(mapindex);
		    for (unsigned int k = 0; k<27; k++) {
		      newmapindex = inputNIt.GetIndex(k);
		      if (newmapindex != mapindex &&
		          input->GetPixel(newmapindex) < dist_act && 
		          input->GetPixel(newmapindex) > -1) {	
			    mapindex_min = newmapindex; 
			    dist_act = input->GetPixel(newmapindex);
		      }	
	        }
			mapindex = mapindex_min;
			new_theta_error = 0;
			new_fi_error = 0;
	      } else {
		    mapindex = newmapindex;
	      }
          // fi:
	      angle_fi_error = new_fi_error;
	      // theta:
	      angle_theta_error = new_theta_error;
		  node_new.SetIndex(mapindex);
		  node_new.SetValue(node.GetValue());
		  // Hay que almacenar el punto real:
		  rn[0]=(float_index[0])*spacing[0]+origin[0];
		  rn[1]=(float_index[1])*spacing[1]+origin[1];
		  rn[2]=(float_index[2])*spacing[2]+origin[2];
		  PointsOutElementFloat.push_back(rn);
		  PointsOutElement.push_back(mapindex);
		  output->SetPixel(mapindex,node.GetValue());
          dist_act = input->GetPixel(mapindex);
          count++;
          if (count > 2000) break;
        } // end while dist_act
		m_PointsOut.push_back(PointsOutElement);
		m_PointsOutFloat.push_back(PointsOutElementFloat);
      } // end if pointsIter 
	} // end if m_seeds 
  } // end if Metodo = 2  
 
  if (m_Metodo == 3) {
	// Implementacion del metodo de Heun 
    float h = m_H;
	float index_temp[3],mapindex_f[3],newmapindex_f[3];
	InputIndexType index_intermedio;
	CovPixelType g_temp,g_aux;
	typedef typename itk::ImageFileWriter< InputType > FileWriter; 
	typename FileWriter::Pointer writer = FileWriter::New();
	
	typename GradRecursiveGaussianFilterType::Pointer grad = GradRecursiveGaussianFilterType::New(); 
	grad->SetInput(input);
	grad->SetSigma(m_Sigma);
    grad->Update();
    
	std::vector<InputIndexType> neigh;
	InputIndexType aux_neigh;
	InputSizeType Size = input->GetLargestPossibleRegion().GetSize();
	if ( m_seeds )  {
  	  typename NodeContainer::ConstIterator pointsIter = m_seeds->Begin();
	  typename NodeContainer::ConstIterator pointsEnd = m_seeds->End();
	  for ( ; pointsIter != pointsEnd; ++pointsIter ) {
	    PointsOutElement.clear();
		PointsOutElementFloat.clear();
	    node = pointsIter.Value();
	    mapindex = node.GetIndex();
        dist_act = input->GetPixel(mapindex);
	    count = 0;
		output->SetPixel(mapindex,2);
	    //std::cout << "Metodo 3. Incluida semilla en " << node.GetIndex() << std::endl;
		mapindex_f[0] = mapindex[0];
		mapindex_f[1] = mapindex[1];
		mapindex_f[2] = mapindex[2];
		int flag_break = 0;
        while (dist_act > 0) {
		  //std::cout << "dist_act " << dist_act << std::endl;
		  //std::cout << "mapindex " << mapindex << std::endl;
		  itk::Vector<float, 3>  G,Gl,Gh;
		  
		  //G[0] = gx->GetOutput()->GetPixel(mapindex);
		  //G[1] = gy->GetOutput()->GetPixel(mapindex);
		  //G[2] = gz->GetOutput()->GetPixel(mapindex);
		  G[0] = grad->GetOutput()->GetPixel(mapindex)[0];
		  G[1] = grad->GetOutput()->GetPixel(mapindex)[1];
		  G[2] = grad->GetOutput()->GetPixel(mapindex)[2];
		  G = G/G.GetNorm();
		  //std::cout << "\nmapindex " << mapindex[0] <<  " " << mapindex[1] << " " << mapindex[2] << std::endl;
		  //std::cout << "mapindex_f " << mapindex_f[0] <<  " " << mapindex_f[1] << " " << mapindex_f[2] << std::endl;
		  index_temp[0] = mapindex_f[0] - h*G[0];
		  index_temp[1] = mapindex_f[1] - h*G[1];
		  index_temp[2] = mapindex_f[2] - h*G[2];
		  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  // Idea original: pto mas cercano
		  if (m_MetodoPtoIntermedio == 1) {
		    index_intermedio[0] = round(index_temp[0]);
            index_intermedio[1] = round(index_temp[1]);
            index_intermedio[2] = round(index_temp[2]);
		    
		    g_temp = grad->GetOutput()->GetPixel(index_intermedio);
		    g_temp[0] = g_temp[0]/g_temp.GetNorm();
		    g_temp[1] = g_temp[1]/g_temp.GetNorm();
		    g_temp[2] = g_temp[2]/g_temp.GetNorm();

		    std::cout << "g_temp " << g_temp[0] << " " << g_temp[1] << " " << g_temp[2] << std::endl;
          }

		  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
		  // Segunda forma para calcular el gradiente en el punto intermedio:: 
		  if (m_MetodoPtoIntermedio == 2) {
		    index_intermedio[0] = round(index_temp[0]);
            index_intermedio[1] = round(index_temp[1]);
            index_intermedio[2] = round(index_temp[2]);
		    float d1,d2;
		  
		    Gl[0] = grad->GetOutput()->GetPixel(index_intermedio)[0];
		    Gl[1] = grad->GetOutput()->GetPixel(index_intermedio)[1];
		    Gl[2] = grad->GetOutput()->GetPixel(index_intermedio)[2];

			Gl = Gl/Gl.GetNorm();
		    d2 = sqrt((index_temp[0] - mapindex_f[0])*(index_temp[0] - mapindex_f[0]) + 
				  (index_temp[1] - mapindex_f[1])*(index_temp[1] - mapindex_f[1]) +
				  (index_temp[2] - mapindex_f[2])*(index_temp[2] - mapindex_f[2]));
			d1 = sqrt((index_temp[0] - index_intermedio[0])*(index_temp[0] - index_intermedio[0]) +
				  (index_temp[1] - index_intermedio[1])*(index_temp[1] - index_intermedio[1]) + 
				  (index_temp[2] - index_intermedio[2])*(index_temp[2] - index_intermedio[2]));
		    g_temp[0] = Gl[0]*d2/(d1+d2) + G[0]*d1/(d1+d2);
		    g_temp[1] = Gl[1]*d2/(d1+d2) + G[1]*d1/(d1+d2);
		    g_temp[2] = Gl[2]*d2/(d1+d2) + G[2]*d1/(d1+d2);
		  
		    //std::cout << "G(mapindex)         " << G << " d2/() " << d2/(d1+d2) << std::endl;
		    //std::cout << "index_intermedio " << index_intermedio[0] <<  " " << index_intermedio[1] << " " << index_intermedio[2] << std::endl;
		    //std::cout << "G(index_intermedio) " << Gl[0] << " " << Gl[1] << " " << Gl[2] << " d1/() " << d1/(d1+d2) << std::endl;
		    //std::cout << "input:  (index_intermedio) " << input->GetPixel(index_intermedio) << " (mapindex)  " << input->GetPixel(mapindex) << std::endl;
		    std::cout << "g_temp " << g_temp[0] << " " << g_temp[1] << " " << g_temp[2] << std::endl;
          }
		  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
		  // Tercera forma, promediando en los vecinos: 
		  if (m_MetodoPtoIntermedio == 3) {
		  neigh.clear();
		  aux_neigh[0]=round(index_temp[0]);
		  aux_neigh[1]=round(index_temp[1]);
		  aux_neigh[2]=round(index_temp[2]);
		  neigh.push_back(aux_neigh);
		  
		  aux_neigh[0]=round(index_temp[0]) + 1;
		  aux_neigh[1]=round(index_temp[1]);
		  aux_neigh[2]=round(index_temp[2]);
		  neigh.push_back(aux_neigh);

          aux_neigh[0]=round(index_temp[0]) - 1;
		  aux_neigh[1]=round(index_temp[1]);
		  aux_neigh[2]=round(index_temp[2]);
		  neigh.push_back(aux_neigh);
          
		  aux_neigh[0]=round(index_temp[0]);
		  aux_neigh[1]=round(index_temp[1]) + 1;
		  aux_neigh[2]=round(index_temp[2]);
		  neigh.push_back(aux_neigh);
          
		  aux_neigh[0]=round(index_temp[0]);
		  aux_neigh[1]=round(index_temp[1]) - 1;
		  aux_neigh[2]=round(index_temp[2]);
          neigh.push_back(aux_neigh); 

          aux_neigh[0]=round(index_temp[0]);
		  aux_neigh[1]=round(index_temp[1]);
		  aux_neigh[2]=round(index_temp[2]) + 1;
          neigh.push_back(aux_neigh); 

          aux_neigh[0]=round(index_temp[0]);
		  aux_neigh[1]=round(index_temp[1]);
		  aux_neigh[2]=round(index_temp[2]) - 1;
          neigh.push_back(aux_neigh); 
		  
		  // Promediamos el valor de gradiente intermedio
		  g_temp.Fill(0);
		  for(unsigned n=0; n<7; n++) {
			if(neigh[n][0]>=0 && neigh[n][1]>=0 && neigh[n][2]>=0 && neigh[n][0]<Size[0] && neigh[n][1]<Size[1] && neigh[n][2]<Size[2]){ 
				float dist=sqrt((index_temp[0]-neigh[n][0])*(index_temp[0]-neigh[n][0])+(index_temp[1]-neigh[n][1])*(index_temp[1]-neigh[n][1])+(index_temp[2]-neigh[n][2])*(index_temp[2]-neigh[n][2]));
				if(input->GetPixel(neigh[n]) != 255 ){
					g_aux = grad->GetOutput()->GetPixel(neigh[n]);
					g_temp = g_temp + (1/(dist+1))*g_aux;
				}
				std::cout << "n = " << n << " dist =  " << dist << " g_temp " <<  g_temp << std::endl; 
			 }
		  }
		  
		  g_temp[0] = g_temp[0]/g_temp.GetNorm();
		  g_temp[1] = g_temp[1]/g_temp.GetNorm();
		  g_temp[2] = g_temp[2]/g_temp.GetNorm();
		           
		  //std::cout << "g_temp " << g_temp[0] << " " << g_temp[1] << " " << g_temp[2] << std::endl;
          
		} // endif metodo_pto_intermedio
		
		  newmapindex_f[0] = mapindex_f[0] - (h/2)*(G[0] + g_temp[0]);
		  newmapindex_f[1] = mapindex_f[1] - (h/2)*(G[1] + g_temp[1]);
		  newmapindex_f[2] = mapindex_f[2] - (h/2)*(G[2] + g_temp[2]);
		  //std::cout << "newmapindex_f " << newmapindex_f[0] <<  " " << newmapindex_f[1] << " " << newmapindex_f[2] << std::endl;
	
		  // Codigo control 
		  float distaux = sqrt((mapindex_f[0]-newmapindex_f[0])*(mapindex_f[0]-newmapindex_f[0]) + 
						   (mapindex_f[1]-newmapindex_f[1])*(mapindex_f[1]-newmapindex_f[1]) +
						   (mapindex_f[2]-newmapindex_f[2])*(mapindex_f[2]-newmapindex_f[2]) );
		  // Codigo control

		  newmapindex[0] = round(newmapindex_f[0]);
		  newmapindex[1] = round(newmapindex_f[1]);
		  newmapindex[2] = round(newmapindex_f[2]);
		  
		  //std::cout << "indextemp " << index_temp[0] << " " << index_temp[1] << " " << index_temp[2] << std::endl;
		  //std::cout << "mapindex " << mapindex <<  " newmapindex " << newmapindex << std::endl;
		  //std::cout << "mapindex_f " << mapindex_f[0] << " " << mapindex_f[1] << " " << mapindex_f[2] <<  " newmapindex_f " << newmapindex_f[0] << " " << newmapindex_f[1] << " " << newmapindex_f[2] << std::endl;
		 
		  if ( input->GetPixel(newmapindex) < 0 ) {
		    //printf("Sale fuera de dominio, decidiendo segun minimo local!\n"); 
		    inputNIt.SetLocation(mapindex);
		    for (unsigned int k = 0; k<27; k++) {
		      newmapindex = inputNIt.GetIndex(k);
		      if (newmapindex != mapindex &&
		          input->GetPixel(newmapindex) < dist_act && 
		          input->GetPixel(newmapindex) > -1) {	
			    mapindex_min = newmapindex; 
			    dist_act = input->GetPixel(newmapindex);
		      }	
	        }
			mapindex = mapindex_min;
			mapindex_f[0] = mapindex_min[0];
		    mapindex_f[1] = mapindex_min[1];
		    mapindex_f[2] = mapindex_min[2];
	      } else {
		    if (mapindex_f == newmapindex_f) {
			  flag_break = 1;
			  std::cout << "mapindex == newmapindex" << std::endl;
			} else {
		      mapindex = newmapindex;
			  mapindex_f[0] = newmapindex_f[0];
		      mapindex_f[1] = newmapindex_f[1];
		      mapindex_f[2] = newmapindex_f[2];
			}
	      }
		  if (flag_break == 0) {
		    //std::cout << "dist_act " << dist_act << std::endl << std::endl;
		    node_new.SetIndex(mapindex);
		    node_new.SetValue(node.GetValue());
		    // Hay que pasar el indicie a coordenadas float:
			rn[0]=(mapindex_f[0])*spacing[0]+origin[0];
			rn[1]=(mapindex_f[1])*spacing[1]+origin[1];
			rn[2]=(mapindex_f[2])*spacing[2]+origin[2];
		    PointsOutElementFloat.push_back(rn);
		
			
		    PointsOutElement.push_back(mapindex);
		    output->SetPixel(mapindex,node.GetValue());
            dist_act = input->GetPixel(mapindex);
			std::cout << "rn " << rn << " distaux " << distaux << " dist_act " << dist_act << std::endl;
            count++;
		  } else {
		    break;
		  }
		  //std::cout << std::endl;
          if (count > 200) break;
		} // end while dist_act
		std::cout << "Num ptos calculados: " << count << std::endl;
		std::cout << "h: " << h << std::endl;
		m_PointsOut.push_back(PointsOutElement);
		m_PointsOutFloat.push_back(PointsOutElementFloat);
	  }
	}
  }
  
  if (m_Metodo == 4) {
	// Implementacion del metodo de Heun (analogo a Runge Kutta de orden)
    float h = m_H;
	float index_temp[3],mapindex_f[3],newmapindex_f[3];
	InputIndexType index_intermedio;
	CovPixelType g_temp,g_aux,K2,K3,K4;
	typedef typename itk::ImageFileWriter< InputType > FileWriter; 
	typename FileWriter::Pointer writer = FileWriter::New();
	
	typename GradRecursiveGaussianFilterType::Pointer grad = GradRecursiveGaussianFilterType::New(); 
	grad->SetInput(input);
	grad->SetSigma(m_Sigma);
    grad->Update();
    
	std::vector<InputIndexType> neigh;
	InputIndexType aux_neigh;
	InputSizeType Size = input->GetLargestPossibleRegion().GetSize();
	if ( m_seeds )  {
  	  typename NodeContainer::ConstIterator pointsIter = m_seeds->Begin();
	  typename NodeContainer::ConstIterator pointsEnd = m_seeds->End();
	  for ( ; pointsIter != pointsEnd; ++pointsIter ) {
	    PointsOutElement.clear();
		PointsOutElementFloat.clear();
	    node = pointsIter.Value();
	    mapindex = node.GetIndex();
        dist_act = input->GetPixel(mapindex);
	    count = 0;
		output->SetPixel(mapindex,2);
	    
		mapindex_f[0] = mapindex[0];
		mapindex_f[1] = mapindex[1];
		mapindex_f[2] = mapindex[2];
		int flag_break = 0;
        while (dist_act > 0) {
		  //std::cout << "dist_act " << dist_act << std::endl;
		  //std::cout << "mapindex " << mapindex << std::endl;
		  itk::Vector<float, 3>  G,Gl,Gh;
		  
		  //G[0] = gx->GetOutput()->GetPixel(mapindex);
		  //G[1] = gy->GetOutput()->GetPixel(mapindex);
		  //G[2] = gz->GetOutput()->GetPixel(mapindex);
		  G[0] = grad->GetOutput()->GetPixel(mapindex)[0];
		  G[1] = grad->GetOutput()->GetPixel(mapindex)[1];
		  G[2] = grad->GetOutput()->GetPixel(mapindex)[2];
		  G = G/G.GetNorm();
		  //std::cout << "\nmapindex " << mapindex[0] <<  " " << mapindex[1] << " " << mapindex[2] << std::endl;
		  //std::cout << "mapindex_f " << mapindex_f[0] <<  " " << mapindex_f[1] << " " << mapindex_f[2] << std::endl;
		  std::cout << "G norm " << G[0] <<  " " << G[1] << " " << G[2] << std::endl;
		  index_temp[0] = mapindex_f[0] - (h/2)*G[0];
		  index_temp[1] = mapindex_f[1] - (h/2)*G[1];
		  index_temp[2] = mapindex_f[2] - (h/2)*G[2];
		  
		  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
		  // Segunda forma para calcular el gradiente en el punto intermedio:: 
		  if (m_MetodoPtoIntermedio >= 1) {
		    index_intermedio[0] = round(index_temp[0]);
            index_intermedio[1] = round(index_temp[1]);
            index_intermedio[2] = round(index_temp[2]);
		    
			K2[0] = grad->GetOutput()->GetPixel(index_intermedio)[0];
		    K2[1] = grad->GetOutput()->GetPixel(index_intermedio)[1];
		    K2[2] = grad->GetOutput()->GetPixel(index_intermedio)[2];
			K2 = K2/K2.GetNorm();
		   
            index_intermedio[0] = round(mapindex_f[0] - (h/2)*K2[0]);
			index_intermedio[1] = round(mapindex_f[1] - (h/2)*K2[1]);
            index_intermedio[2] = round(mapindex_f[2] - (h/2)*K2[2]);
            
			K3[0] = grad->GetOutput()->GetPixel(index_intermedio)[0];
		    K3[1] = grad->GetOutput()->GetPixel(index_intermedio)[1];
		    K3[2] = grad->GetOutput()->GetPixel(index_intermedio)[2];
			K3 = K3/K3.GetNorm();

            /////////
            index_intermedio[0] = round(mapindex_f[0] - (h)*K3[0]);
			index_intermedio[1] = round(mapindex_f[1] - (h)*K3[1]);
            index_intermedio[2] = round(mapindex_f[2] - (h)*K3[2]);
            
			K4[0] = grad->GetOutput()->GetPixel(index_intermedio)[0];
		    K4[1] = grad->GetOutput()->GetPixel(index_intermedio)[1];
		    K4[2] = grad->GetOutput()->GetPixel(index_intermedio)[2];
			K4 = K4/K4.GetNorm();
		    
		  } // endif metodo_pto_intermedio
		  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
		  newmapindex_f[0] = mapindex_f[0] - (h/6)*(G[0] + 2*K2[0] + 2*K3[0] + K4[0]);
		  newmapindex_f[1] = mapindex_f[1] - (h/6)*(G[1] + 2*K2[1] + 2*K3[1] + K4[1]);
		  newmapindex_f[2] = mapindex_f[2] - (h/6)*(G[2] + 2*K2[2] + 2*K3[2] + K4[2]);

		  //std::cout << "newmapindex_f " << newmapindex_f[0] <<  " " << newmapindex_f[1] << " " << newmapindex_f[2] << std::endl;
		  
		  newmapindex[0] = round(newmapindex_f[0]);
		  newmapindex[1] = round(newmapindex_f[1]);
		  newmapindex[2] = round(newmapindex_f[2]);
		  
		  //std::cout << "indextemp " << index_temp[0] << " " << index_temp[1] << " " << index_temp[2] << std::endl;
		  //std::cout << "mapindex " << mapindex <<  " newmapindex " << newmapindex << std::endl;
		  //std::cout << "mapindex_f " << mapindex_f[0] << " " << mapindex_f[1] << " " << mapindex_f[2] <<  " newmapindex_f " << newmapindex_f[0] << " " << newmapindex_f[1] << " " << newmapindex_f[2] << std::endl;
		 
		  if ( input->GetPixel(newmapindex) < 0 ) {
		    //printf("Sale fuera de dominio, decidiendo segun minimo local!\n"); 
		    inputNIt.SetLocation(mapindex);
		    for (unsigned int k = 0; k<27; k++) {
		      newmapindex = inputNIt.GetIndex(k);
		      if (newmapindex != mapindex &&
		          input->GetPixel(newmapindex) < dist_act && 
		          input->GetPixel(newmapindex) > -1) {	
			    mapindex_min = newmapindex; 
			    dist_act = input->GetPixel(newmapindex);
		      }	
	        }
			mapindex = mapindex_min;
			mapindex_f[0] = mapindex_min[0];
		    mapindex_f[1] = mapindex_min[1];
		    mapindex_f[2] = mapindex_min[2];
	      } else {
		    if (mapindex_f == newmapindex_f) {
			  flag_break = 1;
			  std::cout << "mapindex == newmapindex" << std::endl;
			} else {
		      mapindex = newmapindex;
			  mapindex_f[0] = newmapindex_f[0];
		      mapindex_f[1] = newmapindex_f[1];
		      mapindex_f[2] = newmapindex_f[2];
			}
	      }
		  if (flag_break == 0) {
		    //std::cout << "dist_act " << dist_act << std::endl << std::endl;
		    node_new.SetIndex(mapindex);
		    node_new.SetValue(node.GetValue());
		    // Hay que pasar el indicie a coordenadas float:
			rn[0]=(mapindex_f[0])*spacing[0]+origin[0];
			rn[1]=(mapindex_f[1])*spacing[1]+origin[1];
			rn[2]=(mapindex_f[2])*spacing[2]+origin[2];
		    PointsOutElementFloat.push_back(rn);
			std::cout << "rn: " << rn << std::endl; 
		    PointsOutElement.push_back(mapindex);
		    output->SetPixel(mapindex,node.GetValue());
            dist_act = input->GetPixel(mapindex);
            count++;
		  } else {
		    break;
		  }
		  //std::cout << std::endl;
          if (count > 200) break;
		} // end while dist_act
		std::cout << "Num ptos calculados: " << count << std::endl;
		std::cout << "h: " << h << std::endl;
		m_PointsOut.push_back(PointsOutElement);
		m_PointsOutFloat.push_back(PointsOutElementFloat);
	  }
	}
  } // End if metodo == 4
  
}


 // template <class TInput, class TOutput  >
//  VectorContainerType geodesicPath3D<TInput,TOutput>
//  ::GetPointsOut()
//  {
//    return m_PointsOut;
//    
//  } 

  template <class TInput, class TOutput  >
  void geodesicPath3D<TInput,TOutput>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os,indent);
    
  } 



}// end namespaceitk

#endif
