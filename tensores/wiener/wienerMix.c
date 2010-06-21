#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "wienerMix.h"

int main(int argc, char *argv[]) {
  float lambda, gamma;
  unsigned iterations, flag_noise, flag_bias;
  struct vtkFileHeader headerIn,headerOut;
  sscanf(argv[1],"%g",&lambda);
  sscanf(argv[2],"%g", &gamma);
  sscanf(argv[3],"%u",&iterations);
  sscanf(argv[6],"%u",&flag_noise);
  sscanf(argv[7],"%u",&flag_bias);
  strcpy(headerIn.fileName,argv[4]);
  strcpy(headerOut.fileName,argv[5]);
  vtkFileHeaderReader(&headerIn);
  vtkFileDataReader(&headerIn);
  filterData(headerIn,&headerOut,lambda,gamma,iterations,flag_noise, flag_bias);
  vtkFileHeaderWriter(headerOut);
  vtkFileDataWriter(headerOut);
  free(headerIn.data);
  free(headerOut.data);
}

void vtkFileHeaderReader(struct vtkFileHeader *header) {
	FILE *file;
	char line[100],aux[100];
    file=fopen(header->fileName,"r");
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	sscanf(line,"%s %u %u %u",aux,header->dimensions,header->dimensions+1,header->dimensions+2);
	header->size=header->dimensions[0]*header->dimensions[1]*header->dimensions[2];
	fgets(line,1000,file);
	sscanf(line,"%s %g %g %g",aux,header->spacing,header->spacing+1,header->spacing+2);
	fgets(line,1000,file);
	sscanf(line,"%s %g %g %g",aux,header->origin,header->origin+1,header->origin+2);
	fclose(file);
}

void vtkFileDataReader(struct vtkFileHeader *header) {
	FILE *file;
	
	char line[100];
	unsigned i;
	file=fopen(header->fileName,"r");
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	fgets(line,1000,file);
	header->data=(float *)malloc(sizeof(float)*header->size);
	fread((char *)header->data,1,sizeof(float)*header->size,file);

	for(i=0; i<header->size; i++) {
	  swap_4(header->data+i);
	}

	fclose(file);
}

void filterData(struct vtkFileHeader headerIn,struct vtkFileHeader *headerOut,float lambda,float gamma, unsigned iterations, unsigned flag_noise, unsigned flag_bias){

  unsigned count1,count2,count3,count7,count8,count9;
  int aux1,aux2; 
  int count4,count5,count6, countThreshold;
  float *auxIn,*auxMean,*auxVariance,*auxOut;
  float neigh[27],*mean,minVariance,aveVariance,noiseVariance,*copy;

  float edges[6][17]; //six neighborhood edges. 
  float minLocalVariance, meanLocalVariance, localMean[6], localVariance[6];
  unsigned posMinLocalVariance;

  unsigned remainder, quotient;
 
  unsigned ind_edges[6][17]={{9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26},
			     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,14,15,16,17},
			     {1, 2, 4, 5, 7, 8,10,11,14,16,17,19,20,22,23,25,26},
			     {0, 1, 3, 4, 6, 7, 9,10,12,15,16,18,19,21,22,24,25},
			     {0, 1, 2, 3, 4, 5, 9,10,11,12,14,18,19,20,21,22,23},
			     {3, 4, 5, 6, 7, 8,12,14,15,16,17,21,22,23,24,25,26}};
 
  FILE *file;
  float SNR, localB, localMS, x_0, x_1, *localS, *SNR_table;  
  
 
  clock_t time, timeAux;

  // headerOut
  headerOut->dimensions[0]=headerIn.dimensions[0];
  headerOut->dimensions[1]=headerIn.dimensions[1];
  headerOut->dimensions[2]=headerIn.dimensions[2];
  headerOut->spacing[0]=headerIn.spacing[0];
  headerOut->spacing[1]=headerIn.spacing[1];
  headerOut->spacing[2]=headerIn.spacing[2];
  headerOut->origin[0]=headerIn.origin[0];
  headerOut->origin[1]=headerIn.origin[1];
  headerOut->origin[2]=headerIn.origin[2];
  headerOut->size=headerIn.size;
  headerOut->data=(float *)malloc(sizeof(float)*headerOut->size);
  mean=(float *)malloc(sizeof(float)*headerOut->size);
  copy=(float *)malloc(sizeof(float)*headerOut->size);
  localS=(float*)malloc(sizeof(float)*headerOut->size);

  for (count1=0;count1<headerOut->size;count1++){
      copy[count1]=headerIn.data[count1];
      headerOut->data[count1]=copy[count1];
  }
  
  aux1=headerIn.dimensions[1]*headerIn.dimensions[0];
  aux2=headerIn.dimensions[0];
  
	
 
 //Remove bias
  if(flag_bias==1){
    SNR_table=(float *)malloc(sizeof(float)*10001);
    file=fopen("SNR.dat","r");
    fread((char *)SNR_table,1,sizeof(float)*10001,file);
    fclose(file);    
    
    auxIn=copy+aux1;
    auxMean=mean+aux1;
    auxVariance=headerOut->data+aux1;

    count8=0;
    for (count3=1;count3<headerOut->dimensions[2]-1;count3++) {
      auxIn+=aux2;
      auxMean+=aux2;
      auxVariance+=aux2;
      for (count2=1;count2<headerOut->dimensions[1]-1;count2++) {
	auxIn++;
	auxMean++;
	auxVariance++;
	for (count1=1;count1<headerOut->dimensions[0]-1;count1++) {
	  count7=0;
	  for (count4=-1;count4<=1;count4++) {
	    for (count5=-1;count5<=1;count5++) {
	      for (count6=-1;count6<=1;count6++) {
		
		neigh[count7]=(*(auxIn+count4+(count5*aux1)+(count6*aux2))); 
		
		count7++;
	      }
	    }
	  }

	  // Edges of neighborhood
	  for(count5=0; count5<6; count5++){
	    for (count4=0; count4<17; count4++){
	      edges[count5][count4]=neigh[ind_edges[count5][count4]];
	    }
	  }
	  
	  *auxMean=0.0;
	  for (count4=0;count4<27;count4++) {
	    *auxMean+=neigh[count4]/27.0f;
	  }
	  
	  *auxVariance=0.0;
	  for (count4=0;count4<27;count4++) {
	    *auxVariance+=(neigh[count4]-*auxMean)*(neigh[count4]-*auxMean)/26.0;
	  }

	  minLocalVariance=0.0;
	  
	  for (count4=0; count4<6; count4++){
	    // Edge Mean
	    localMean[count4]=0.0;
	    localVariance[count4]=0.0;
	    for (count5=0; count5<17; count5++) {
	      localMean[count4]+= edges[count4][count5]/17.0;
	    }

	    // EdgeVariance
	    for (count5=0; count5<17; count5++) {
	      localVariance[count4]+=(edges[count4][count5]-localMean[count4])*(edges[count4][count5]-localMean[count4])/16.0;
	    }
	    
	    if (localVariance[count4]<minLocalVariance || minLocalVariance==0){
	      minLocalVariance=localVariance[count4];
	      posMinLocalVariance=count4;
	    }
	  }
	  
	  if(*auxVariance>minLocalVariance){
	    *auxMean=localMean[posMinLocalVariance];
	    *auxVariance=minLocalVariance;
	  }
  
	  *auxVariance=sqrt(*auxVariance);
	  SNR=(*auxMean)/(*auxVariance);
	  
	  if(SNR<=SNR_table[0]) {
	    localS[count8]=0.0;
	  }
	  else{
	    localMS=(*auxMean)*(*auxMean)+(*auxVariance)*(*auxVariance);
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
	    localS[count8]=localB*sqrt(localMS/(2+(localB*localB)));
	  }
	  
	  count8++;
	  auxIn++;
	  auxMean++;
	  auxVariance++;
	}
	auxIn++;
	auxMean++;
	auxVariance++;
      }
      auxIn+=aux2;
      auxMean+=aux2;
      auxVariance+=aux2;
    }
    
    count8=0;
    auxIn=copy+aux1;
    auxMean=mean+aux1;
 
    for (count3=1;count3<headerOut->dimensions[2]-1;count3++) {
      auxIn+=aux2;
      auxMean+=aux2;
      for (count2=1;count2<headerOut->dimensions[1]-1;count2++) {
	auxIn++;
	auxMean++;
	for (count1=1;count1<headerOut->dimensions[0]-1;count1++) {
	  *auxIn=(*auxIn)-(*auxMean)+localS[count8];
	  auxIn++;
	  auxMean++;
	  count8++;

	}
	auxIn++;
	auxMean++;
      }
      auxIn+=aux2;
      auxMean+=aux2;
    }
  }
  //End remove bias

  timeAux=clock();

  //ITERATIONS 
  for (count9=0;count9<iterations;count9++) {
    auxIn=copy+aux1;
    auxMean=mean+aux1;
    auxVariance=headerOut->data+aux1;

    minVariance=0.0;
    aveVariance=0.0;
	
    count8=0;
   
    // for Z
    for (count3=1;count3<headerOut->dimensions[2]-1;count3++) {
      auxIn+=aux2;
      auxMean+=aux2;
      auxVariance+=aux2;

      // for Y
      for (count2=1;count2<headerOut->dimensions[1]-1;count2++) {
	auxIn++;
	auxMean++;
	auxVariance++;

	// for X
	for (count1=1;count1<headerOut->dimensions[0]-1;count1++) {
	  count7=0;

	  // Build the neighborhood
	  for (count4=-1;count4<=1;count4++) {
	    for (count5=-1;count5<=1;count5++) {
	      for (count6=-1;count6<=1;count6++) {
       		neigh[count7]=(*(auxIn+count4+(count5*aux1)+(count6*aux2))); 
		count7++;

	      }
	    }
	  }//end build neighborhood

	  // Edges of neighborhood
	  for(count5=0; count5<6; count5++){
	    for (count4=0; count4<17; count4++){
	      edges[count5][count4]=neigh[ind_edges[count5][count4]];
	    }
	  }
	 
	  // End of building the six-subneighborhoods
	  
	  *auxMean=0.0;

	  minLocalVariance=0.0;

	  // Mean in 3x3x3 neighboorhood
	  for (count4=0; count4<27; count4++){
	    *auxMean+=neigh[count4]/27.0;
	  }
	  *auxVariance=0.0;
	  // Variance in 3x3x3 neighborhood
	  for (count4=0; count4<27; count4++){
	    *auxVariance+=(neigh[count4]-*auxMean)*(neigh[count4]-*auxMean)/26.0;
	  }
	
	  //localMean and localVariance in the 6 edges neighborhood
	  for (count4=0; count4<6; count4++){
	    
	    localMean[count4]=0.0;
	    localVariance[count4]=0.0;
	    
	    // 6 local means
	    for (count5=0; count5<17; count5++) {
	      localMean[count4]+= edges[count4][count5]/17.0;
	    }

	    // 6 local variances
	    for (count5=0; count5<17; count5++) {
	      localVariance[count4]+=(edges[count4][count5]-localMean[count4])*(edges[count4][count5]-localMean[count4])/16.0;
	    }

	    //minLocalVariance: minimal variance between the 6 edge neighborhood
	    if (localVariance[count4]<minLocalVariance || minLocalVariance==0){
	      minLocalVariance=localVariance[count4];
	      posMinLocalVariance=count4;  //MinVarianceIndexes
	    }
	  }
	  
	  
	  //MinLocalMean=localMean[posMinLocalVariance];
	  
	  meanLocalVariance=0.0;
	  //meanLocalVariance: mean of the variance in the 6 edges neighborhood.
	  for(count4=0; count4<6; count4++)
	  {
	    meanLocalVariance+=localVariance[count4]/6.0;
	  }

	  //ThresoldedIndexes
	  if(((meanLocalVariance-minLocalVariance)/(meanLocalVariance))>gamma){
	    *auxMean=localMean[posMinLocalVariance];  
	    *auxVariance=minLocalVariance;
	  }
	  
  
	  if (minVariance>*auxVariance || minVariance == 0) {
	    minVariance=*auxVariance;
	  }
	  aveVariance+=*auxVariance;
	  count8++;
	  auxIn++;
	  auxMean++;
	  auxVariance++;

	}// X
	auxIn++;
	auxMean++;
	auxVariance++;
      } //Y
      auxIn+=aux2;
      auxMean+=aux2;
      auxVariance+=aux2;
    }
    

    aveVariance/=((float)count8);
    noiseVariance=minVariance*(1-lambda)+aveVariance*lambda;
    printf("%g\n", noiseVariance);
  
    auxIn=copy+aux1;
    auxMean=mean+aux1;
    auxVariance=headerOut->data+aux1;
    auxOut=headerOut->data+aux1;
    for (count3=1;count3<headerOut->dimensions[2]-1;count3++) {
      auxIn+=aux2;
      auxMean+=aux2;
      auxVariance+=aux2;
      auxOut+=aux2;
      for (count2=1;count2<headerOut->dimensions[1]-1;count2++) {
	auxIn++;
	auxMean++;
	auxVariance++;
	auxOut++;
	for (count1=1;count1<headerOut->dimensions[0]-1;count1++) {
	  if(flag_noise==1){
	    *auxVariance=*auxVariance-noiseVariance;
	    if(*auxVariance<0){
	      *auxVariance=0.0;
	    }
	  }

	  *auxOut=(*auxIn-*auxMean)*(*auxVariance)/(*auxVariance+noiseVariance)+(*auxMean);
	  auxIn++;
	  auxMean++;
	  auxVariance++;
	  auxOut++;
	}
	auxIn++;
	auxMean++;
	auxVariance++;
	auxOut++;
      }
      auxIn+=aux2;
      auxMean+=aux2;
      auxVariance+=aux2;
      auxOut+=aux2;
    }
    for (count1=0;count1<headerOut->size;count1++) {
      copy[count1]=headerOut->data[count1];
    }
    time=clock();
    printf("ITERATION %d: TIME %f\n",count9+1,(time-timeAux)/((float) CLOCKS_PER_SEC));
    timeAux=time;
  }
  free(mean);
  free(copy);
}

void vtkFileHeaderWriter(struct vtkFileHeader header) {
	FILE *file;
    file = fopen(header.fileName,"w");
	fprintf(file,"# vtk DataFile Version 3.0\n");
	fprintf(file,"vtk output\n");
	fprintf(file,"BINARY\n");
	fprintf(file,"DATASET STRUCTURED_POINTS\n");
	fprintf(file,"DIMENSIONS %u %u %u\n",header.dimensions[0],header.dimensions[1],header.dimensions[2]);
	fprintf(file,"SPACING %f %f %f\n",header.spacing[0],header.spacing[1],header.spacing[2]);
	fprintf(file,"ORIGIN %f %f %f\n",header.origin[0],header.origin[1],header.origin[2]);
	fprintf(file,"POINT_DATA %u\n",header.size);
	fprintf(file,"SCALARS scalars float 1\n");
	fprintf(file,"LOOKUP_TABLE default\n");
	fclose(file);
}

void vtkFileDataWriter(struct vtkFileHeader header) {
	FILE *file;
	unsigned i;

	for(i=0; i<header.size; i++) {
	  swap_4(&header.data[i]);
	}
	file = fopen(header.fileName,"a");
	fwrite((char *) header.data,1,sizeof(float)*header.size,file);
	fclose(file);
}

void swap_4(void *tnf4)              /* 4 byte floating point numbers */
{
 int *tni4=(int *)tnf4;
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}
