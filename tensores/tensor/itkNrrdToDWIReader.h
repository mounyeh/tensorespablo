/*=========================================================================

  Program:   itkNrrdToDWIReader.h
  Module:    $$
  Language:  C++
  Date:      $Date: 2008/06/04  $
  Version:   $ $

=========================================================================*/
#ifndef __itkNrrdToDWIReader_h
#define __itkNrrdToDWIReader_h

// Undefine an eventual DWImages macro
#ifdef NrrdToDWIReader
#undef NrrdToDWIReader
#endif
#include <itkImage.h>
#include "itkDWImages.h"
#include <itkImageFileReader.h>


namespace itk
{

template < class TOutput >
class NrrdToDWIReader : public itk::ImageFileReader< TOutput >
{
public:
	/** Standard class typedefs. */
	typedef NrrdToDWIReader						    Self;
	typedef itk::ImageFileReader<TOutput>			Superclass;
	typedef SmartPointer<Self>						Pointer;
	
	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro( NrrdToDWIReader, ImageFileReader );
	
	/** Propagating some typedef from the superclass */
	typedef TOutput										OutputImageType;
	typedef typename TOutput::SizeType					SizeType;
	typedef typename TOutput::SpacingType				SpacingType;
	typedef typename TOutput::DirectionType				DirectionType;
	typedef typename TOutput::RegionType				ImageRegionType;
	typedef typename TOutput::IndexType					IndexType;
	typedef typename TOutput::DirectionVectorType       DirectionVectorType;
	typedef typename TOutput::DiffusionDirectionsType	DiffusionDirectionsType;
	typedef typename TOutput::IndicatorType             IndicatorType;
	typedef typename TOutput::BValuesType               BValuesType;

	void GenerateData();
protected:
	/** Default Constructor. */
	NrrdToDWIReader(){}
	~NrrdToDWIReader(){}
};

} // end namespace itk


#if ITK_TEMPLATE_TXX
# include "itkNrrdToDWIReader.txx"
#endif

#endif
