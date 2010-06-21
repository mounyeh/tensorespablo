/*=========================================================================

  Program:   itkDICOMtoDWIReader.h
  Module:    $$
  Language:  C++
  Date:      $Date: 2008/06/04  $
  Version:   $ $

=========================================================================*/
#ifndef __itkDICOMtoDWIReader_h
#define __itkDICOMtoDWIReader_h

// Undefine an eventual DWImages macro
#ifdef DICOMtoDWIReader
#undef DICOMtoDWIReader
#endif
#include <itkImage.h>
#include <itkImageSource.h>
#include "itkGDCMSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkArray.h"


namespace itk
{

template < class TOutput >
class DICOMtoDWIReader : public itk::ImageSource<TOutput>
{
public:
	/** Standard class typedefs. */
	typedef DICOMtoDWIReader						Self;
	typedef itk::ImageSource<TOutput>				Superclass;
	typedef SmartPointer<Self>						Pointer;
	
	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(ImageFileReader, ImageSource);
	
	/** Specify the file to read. This is forwarded to the IO instance. */
	itkSetStringMacro(DirectoryName);
	itkGetStringMacro(DirectoryName);

	/** Propagating some typedef from the superclass */
	
	typedef TOutput											TOutputType;
	typedef typename TOutputType::SizeType					SizeType;
	typedef typename TOutputType::SpacingType				SpacingType;
	typedef typename TOutputType::DirectionType				DirectionType;
	typedef typename TOutputType::RegionType				VectorImageRegionType;
	typedef typename TOutputType::DiffusionDirectionsType	DiffusionDirectionsType;
	typedef typename TOutputType::IndexType					IndexType;

	typedef itk::Image<float, 3>						ImageType;

	typedef itk::GDCMSeriesFileNames					InputNamesGeneratorType;
	typedef itk::ImageSeriesReader< ImageType>		    ImageReaderType;

	bool SliceOrderIS;

	/** Prepare the allocation of the output image during the first back
	* propagation of the pipeline. */
	virtual void GenerateOutputInformation(void);

	void GenerateData();

	void InsertUnique( std::vector<float> & vec, float value )	
	{
		int n = vec.size();
		if (n == 0)
		{
			vec.push_back( value );
			return;
		}

		for (int k = 0; k < n ; k++)
		{
			if (vec[k] == value)
			{
				return;
			}
		}

		// if we get here, it means value is not in vec.
		vec.push_back( value );
		return;

	};
	 
protected:
	/** Default Constructor. */
	DICOMtoDWIReader();
	~DICOMtoDWIReader();
	
	std::string  m_DirectoryName;
};

} // end namespace itk


#if ITK_TEMPLATE_TXX
# include "itkDICOMtoDWIReader.txx"
#endif

#endif
