/*=========================================================================

  Program:   itkDWImages
  Module:    $ $
  Language:  C++
  Date:      $Date: 2008/06/04$
  Version:   $ $

=========================================================================*/
#ifndef _itkDWImages_txx
#define _itkDWImages_txx

#include "itkDWImages.h"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

namespace itk
{
/*---------------------------------------------------------------------------**/
/** Constructors...*/
template<class T, unsigned int D>
DWImages<T,D>::DWImages():Superclass(){
	m_NumImages    = 7;
	m_NumBaselines = 1;
	m_NumDWImages  = 6;	
}

template<class T, unsigned int D>
DWImages<T,D>::~DWImages( )
{
}

template<class T, unsigned int D> 
void DWImages<T,D>::GetImageComponent( unsigned int index, ScalarImagePointerType image){
	
	typedef itk::ImageRegionIteratorWithIndex< Self >  VectorIteratorType;
	typedef itk::ImageRegionIterator< ScalarImageType> ImagesIteratorType;

	VectorIteratorType it1(this, this->GetRequestedRegion());
	ImagesIteratorType it2(image, image->GetRequestedRegion());

	for ( it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1,++it2 )
		it2.Set(   ( it1.Get() )[index]   );
}


template<class T, unsigned int D> 
void DWImages<T,D>::GetT2Image(ScalarImagePointerType image){
	GetImageComponent( m_IndexesBaselines[0], image);
}


template<class T, unsigned int D> 
void DWImages<T,D>::CreateDiffusionDirections(DiffusionDirectionsType directions){
	const float eps = 1e-10;
	DirectionVectorType aux_vector;
	m_NumBaselines = 0;
	m_NumDWImages  = 0;
	
	bool* baselines = new bool[m_NumImages];
	for( unsigned int i=0; i<m_NumImages; ++i){
		baselines[i] = true;
		aux_vector[0]=0.0;
		aux_vector[1]=0.0;
		aux_vector[2]=0.0;
		
		for( unsigned int j=0; j<3; ++j ){
			// If is not a baseline
			//std::cout << "i: " << i << "  j: " << j << "  directions[i][j]: " << directions[i][j] << std::endl;
			if(       !(directions[i][j]!=directions[i][j])   &&   fabs(directions[i][j])>eps   ){
				aux_vector[j] = directions[i][j];
				baselines[i] = false;
			}
		}
		if( baselines[i] )
			m_NumBaselines++;
		else{
			m_NumDWImages++;
			for(unsigned int j=0; j<3; ++j){
				if(aux_vector[j]!=aux_vector[j]){
					aux_vector[j]=0.0;
				}
			}
			m_DiffusionDirections.push_back(aux_vector);
		}
	}
	unsigned int p1=0, p2=0;
	m_IndexesDWImages.SetSize( m_NumDWImages );
	m_IndexesBaselines.SetSize( m_NumBaselines );
	for( unsigned int k=0; k<m_NumImages; ++k ){
		if( baselines[k] )
			m_IndexesBaselines[p1++] = k;
		else
			m_IndexesDWImages[p2++]  = k;
	}
	delete[] baselines;
}



template<class T, unsigned int D>
void DWImages<T,D>::SetDiffusionDirections(DiffusionDirectionsType directions){
	m_DiffusionDirections.clear();
	for( unsigned int k=0; k<directions.size(); ++k )
		m_DiffusionDirections.push_back(directions[k]);
}

template<class T, unsigned int D> 
void DWImages<T,D>::GetDiffusionDirections( DiffusionDirectionsType& directions ){
	directions.clear();
	for(unsigned int i=0; i<m_NumDWImages; i++)
		directions.push_back( m_DiffusionDirections[i] );
}

 	
} // end namespace itk

#endif
