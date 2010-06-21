/*=========================================================================

  Program:  itkDWImages
  Module:    $ $
  Language:  C++
  Date:      $Date: 2008/06/04$
  Version:   $ $

=========================================================================*/
#ifndef __itkDWImages_h
#define __itkDWImages_h

// Undefine an eventual DWImages macro
#ifdef DWImages
#undef DWImages
#endif

#include <itkImage.h>
#include <itkObject.h>
#include <itkVectorImage.h>
#include <itkArray.h>
#include "itkVectorImageNeighborhoodAccessorFunctor.h"
#include "itkDefaultVectorPixelAccessor.h"
#include "itkDefaultVectorPixelAccessorFunctor.h"

namespace itk
{

template < typename TComponent, unsigned int NDimension=3>
class DWImages : public itk::VectorImage < TComponent,  NDimension>
{
public:
	/** Standard class typedefs. */
	typedef		DWImages										Self;
	typedef		itk::VectorImage < TComponent,  NDimension>		Superclass;
	typedef		SmartPointer<Self>								Pointer;
	typedef     SmartPointer<const Self>                        ConstPointer;
	typedef     WeakPointer<const Self>                         ConstWeakPointer;
	
	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro( DWImages, VectorImage );
	/*================================================================================================*/
	/** Neighborhood acces to DWImages*/
	typedef		VectorImageNeighborhoodAccessorFunctor<Self>    NeighborhoodAccessorFunctorType;
	/** Accessor type that convert data between internal and external representations.  */
	typedef DefaultVectorPixelAccessor< TComponent >            AccessorType;
	/** Functor to provide a common API between DefaultPixelAccessor and DefaultVectorPixelAccessor */
	typedef DefaultVectorPixelAccessorFunctor< Self >           AccessorFunctorType;
	/** Return the Pixel Accessor object */
	AccessorType GetPixelAccessor( void ){ return AccessorType( this->GetVectorLength() ); }
	/** Return the Pixel Accesor object */
	const AccessorType GetPixelAccessor( void ) const{ return AccessorType( this->GetVectorLength() ); }
	NeighborhoodAccessorFunctorType	GetNeighborhoodAccessor(){
		return NeighborhoodAccessorFunctorType(  this->GetVectorLength()  ); 
	}
	const NeighborhoodAccessorFunctorType	GetNeighborhoodAccessor() const {
		return NeighborhoodAccessorFunctorType(  this->GetVectorLength()  ); 
	}
	/*================================================================================================*/
	
	typedef		TComponent						ComponentType;
	
	typedef typename itk::Image<ComponentType, NDimension>			ScalarImageType;
	typedef typename ScalarImageType::Pointer				ScalarImagePointerType;
	typedef typename ScalarImageType::RegionType				ScalarImageRegionType;

	typedef	typename itk::Vector<double,3>				        DirectionVectorType;
	typedef			 std::vector< DirectionVectorType>		DiffusionDirectionsType;
	
	typedef itk::Array<unsigned int>                                	IndicatorType;
	typedef itk::Array<double>                                      	BValuesType;
	
	
	
	itkSetMacro(NumImages, unsigned int); 
	itkGetMacro(NumImages, unsigned int); 
	
	itkSetMacro(NumBaselines, unsigned int); 
	itkGetMacro(NumBaselines, unsigned int); 
	
	itkSetMacro(NumDWImages, unsigned int); 
	itkGetMacro(NumDWImages, unsigned int); 
		
	itkSetMacro(B_Value, double); 
	itkGetMacro(B_Value, double); 

	itkSetMacro(UniqueB, bool); 
	itkGetMacro(UniqueB, bool); 

	void SetIndexesBaselines(IndicatorType indexes){
		m_IndexesBaselines = indexes;
	}
	IndicatorType GetIndexesBaselines( void ){
		return m_IndexesBaselines;
	}
	
	void SetIndexesDWImages(IndicatorType indexes){
		m_IndexesDWImages = indexes;
	}
	
	IndicatorType GetIndexesDWImages( void ){
		return m_IndexesDWImages;
	}
	
	void SetBValues(BValuesType values){
		const double tol = 1.0f;
		m_BValues = values;
		m_B_Value  = values[0];
		m_UniqueB = true;
		if( values.GetSize()>1 ){
			double* nonZeroValues = new double[values.GetSize() ];
			nonZeroValues[0]      = m_B_Value;
			unsigned int pos = 0;
			for( unsigned int k=0; k<values.GetSize(); ++k ){
				if( values[k]>0.5f )
					nonZeroValues[pos++] = values[k];
			}
			m_B_Value = nonZeroValues[0];
			for( unsigned int k=1; k<pos; ++k ){
				if(   fabs( nonZeroValues[k-1] - nonZeroValues[k] )   >   tol   ){
					m_UniqueB = false;
					break;
				}
			}
			delete[] nonZeroValues;
		}
	}
	
	BValuesType GetBValues( void ){
		return m_BValues;
	}
	
	void CreateDiffusionDirections(DiffusionDirectionsType);	
	void SetDiffusionDirections(DiffusionDirectionsType);
	void GetDiffusionDirections(DiffusionDirectionsType& directions );
	
	void GetImageComponent(unsigned int, ScalarImagePointerType); 
		
	void GetT2Image(ScalarImagePointerType); 


protected:

	/** Default Constructor. */
	DWImages();
	~DWImages();
	
	unsigned int                m_NumImages;
	unsigned int                m_NumBaselines;
	unsigned int                m_NumDWImages;
	DiffusionDirectionsType     m_DiffusionDirections;
	IndicatorType               m_IndexesDWImages;
	IndicatorType               m_IndexesBaselines;
	BValuesType                 m_BValues;   // BValues vector
	double                      m_B_Value;	 // One BValue vector, if it is the same for every acquisition
	bool                        m_UniqueB;	 // true if the bValue is the same for every DWI
};

} // end namespace itk

#if ITK_TEMPLATE_TXX
# include "itkDWImages.txx"
#endif

#endif
