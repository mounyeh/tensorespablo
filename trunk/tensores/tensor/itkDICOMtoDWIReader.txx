/*=========================================================================

  Program:	 
  Module:    DICOMtoDWIReader.txx
  Language:  C++
  Date:      $Date: 2008/06/04 $
  Version:   $ $

=========================================================================*/
#ifndef _itkDICOMtoDWIReader_txx
#define _itkDICOMtoDWIReader_txx

#include "itkDICOMtoDWIReader.h"
#include <itkImageSeriesReader.h>
#include "gdcmDictSet.h"        // access to dictionary

#include "gdcmDictEntry.h"      // access to dictionary
#include "gdcmGlobal.h"         // access to dictionary


namespace itk
{
/*---------------------------------------------------------------------------**/
/** Constructors...*/
template<class T>
DICOMtoDWIReader<T>::DICOMtoDWIReader():Superclass(){
	m_DirectoryName="";
	const gdcm::DictEntry GEDictBValue( 0x0043, 0x1039, "IS", "1", "B Value of diffusion weighting" );
	const gdcm::DictEntry GEDictXGradient( 0x0019, 0x10bb, "DS", "1", "X component of gradient direction" );	
	const gdcm::DictEntry GEDictYGradient( 0x0019, 0x10bc, "DS", "1", "Y component of gradient direction" );
	const gdcm::DictEntry GEDictZGradient( 0x0019, 0x10bd, "DS", "1", "Z component of gradient direction" );
	const gdcm::DictEntry GEDictZOrientation( 0x0019, 0x101a, "LO", "1", "Orientation of Z axis" );

	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictBValue);
	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictXGradient);
	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictYGradient);
	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictZGradient);
	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictZOrientation);
	
};


template <class T>
DICOMtoDWIReader<T>::~DICOMtoDWIReader()
{
}


template <class T>
void 
DICOMtoDWIReader<T>::GenerateOutputInformation(void)
{
	typename T::Pointer output = this->GetOutput();
	
	itkDebugMacro(<<"Reading directory for GenerateOutputInformation()" << m_DirectoryName);

	// Check if the DirectoryName is ok
	InputNamesGeneratorType::Pointer inputNames = InputNamesGeneratorType::New();
	inputNames->SetInputDirectory( m_DirectoryName );

	ImageReaderType::FileNamesContainer filenames; 
  	
	const ImageReaderType::FileNamesContainer & filenamesInSeries = inputNames->GetInputFileNames();
	  
	// If there are just one file in the series returned above, it is obvious that the series is 
	// not complete. There is no way to put diffusion weighted images in one file, even for mosaic 
	// format. We, then, need to find all files in that directory and treate them as a single series
	// of diffusion weighted image.

	if ( filenamesInSeries.size() > 1 ) 
	{ 
		int nFiles = filenamesInSeries.size(); 
		filenames.resize( 0 ); 
		for (int k = 0; k < nFiles; k++) 
		{ 
			filenames.push_back( filenamesInSeries[k] ); 
		} 
	} 
	else 
	{ 
		std::cout << "gdcm returned just one file. \n"; 
		return;
	}
   
	typedef itk::GDCMImageIO ImageIOType;
	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	 
	ImageReaderType::Pointer reader = ImageReaderType::New();
	reader->SetImageIO( gdcmIO );
	reader->SetFileNames( filenames );
	try
	{
      reader->Update();
    }
	catch (itk::ExceptionObject &excp)
    {
      std::cerr << "Exception thrown while reading the series" << std::endl;
      std::cerr << excp << std::endl;
      return;
    }

	//const gdcm::DictEntry GEDictBValue( 0x0043, 0x1039, "IS", "1", "B Value of diffusion weighting" );
//	const gdcm::DictEntry GEDictXGradient( 0x0019, 0x10bb, "DS", "1", "X component of gradient direction" );	
//	const gdcm::DictEntry GEDictYGradient( 0x0019, 0x10bc, "DS", "1", "Y component of gradient direction" );
//	const gdcm::DictEntry GEDictZGradient( 0x0019, 0x10bd, "DS", "1", "Z component of gradient direction" );
//
//	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictBValue);
//	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictXGradient);
//	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictYGradient);
//	gdcm::Global::GetDicts()->GetDefaultPubDict()->AddEntry(GEDictZGradient);

	ImageReaderType::DictionaryArrayRawPointer inputDict = reader->GetMetaDataDictionaryArray();
	std::string tag;
	
	//-----------------------------------------------------------------------------
	// Image Size
	SizeType dimSize;
	tag.clear();
	itk::ExposeMetaData<std::string> ( *(*inputDict)[0], "0028|0010", tag );
	dimSize[0] = atoi( tag.c_str() );
	
	tag.clear();
	itk::ExposeMetaData<std::string> ( *(*inputDict)[0], "0028|0011", tag );
	dimSize[1] = atoi( tag.c_str() );
	
	int nSlice= inputDict->size();
	
	std::vector<float> sliceLocations(0);
	for (int k = 0; k < nSlice; k++)
    {
      tag.clear();
      itk::ExposeMetaData<std::string> ( *(*inputDict)[k], "0020|1041",  tag);
      float sliceLocation = atof( tag.c_str() );
      InsertUnique( sliceLocations, sliceLocation );
    }    
	
	
	dimSize[2] = sliceLocations.size();
	
	output->SetLargestPossibleRegion(dimSize);
	//--------------------------------------------------------------------
	
	// Image Spacing:
	float fspacing[3];
	SpacingType spacing;
	tag.clear();
	itk::ExposeMetaData<std::string> ( *(*inputDict)[0], "0028|0030", tag );
	sscanf( tag.c_str(), "%f\\%f", &fspacing[0], &fspacing[1] );
	
	
	tag.clear();
	itk::ExposeMetaData<std::string> ( *(*inputDict)[0], "0018|0088", tag );
	fspacing[2] = atof( tag.c_str() );

	spacing[0]=fspacing[0];
	spacing[1]=fspacing[1];
	spacing[2]=fspacing[2];
	output->SetSpacing(spacing);
	//---------------------------------------------------------------------
	
	// Image Origin:
	float origin[3];
	tag.clear();
	itk::ExposeMetaData<std::string> ( *(*inputDict)[0], "0020|0032", tag );
	sscanf( tag.c_str(), "%f\\%f\\%f", &origin[0], &origin[1], &origin[2] );

	output->SetOrigin(origin);
	//---------------------------------------------------------------------
	
	// Adapt to RAS reference system
	std::vector<double> axis;
	
	// check ImageOrientationPatient and figure out slice direction in
	// L-P-I (right-handed) system.
	// In Dicom, the coordinate frame is L-P by default. Look at
	// http://medical.nema.org/dicom/2007/07_03pu.pdf ,  page 301
	tag.clear();
	itk::ExposeMetaData<std::string> ( *(*inputDict)[0], "0020|0037", tag );
	float xRow, yRow, zRow, xCol, yCol, zCol, xSlice, ySlice, zSlice;
	sscanf( tag.c_str(), "%f\\%f\\%f\\%f\\%f\\%f", &xRow, &yRow, &zRow, &xCol, &yCol, &zCol );
	
	// In Dicom, the measurement frame is L-P by default. Look at
	// http://medical.nema.org/dicom/2007/07_03pu.pdf ,  page 301, in
	// order to make this compatible with Slicer's RAS frame, we
	// multiply the direction cosines by the negatives of the resolution
	// (resolution is required by nrrd format). Direction cosine is not
	// affacted since the resulting frame is still a right-handed frame.
	xRow = -xRow;
	yRow = -yRow;

	xCol = -xCol;
	yCol = -yCol;

	// Cross product, this gives I-axis direction
	xSlice = (yRow*zCol-zRow*yCol)*spacing[2];
	ySlice = (zRow*xCol-xRow*zCol)*spacing[2];
	zSlice = (xRow*yCol-yRow*xCol)*spacing[2];
	
	xRow *= spacing[0];
	yRow *= spacing[0];
	zRow *= spacing[0];

	xCol *= spacing[1];
	yCol *= spacing[1];
	zCol *= spacing[1];
	
	float x0, y0, z0;
	float x1, y1, z1;
	tag.clear();
	itk::ExposeMetaData<std::string> ( *(*inputDict)[0], "0020|0032", tag );
	sscanf( tag.c_str(), "%f\\%f\\%f", &x0, &y0, &z0 );
	tag.clear();
	
	SliceOrderIS=true;
	// assume volume interleaving, i.e. the second dicom file stores
	// the second slice in the same volume as the first dicom file
	itk::ExposeMetaData<std::string> ( *(*inputDict)[1], "0020|0032", tag );
	sscanf( tag.c_str(), "%f\\%f\\%f", &x1, &y1, &z1 );
	x1 -= x0; y1 -= y0; z1 -= z0;
	x0 = x1*xSlice + y1*ySlice + z1*zSlice;
	if (x0 < 0)
	{
		SliceOrderIS = false;
	}
	
	DirectionType direction;
	
	direction[0][0]=xRow;
	direction[1][0]=yRow;
	direction[2][0]=zRow;

	direction[0][1]=xCol;
	direction[1][1]=yCol;
	direction[2][1]=zCol;	
	
	if(SliceOrderIS){
		direction[0][2]=xSlice;
		direction[1][2]=ySlice;
		direction[2][2]=zSlice;
	}else{
		
		origin[2]=dimSize[2]*spacing[2]-origin[2];
		output->SetOrigin(origin);
		
		direction[0][2]=-xSlice;
		direction[1][2]=-ySlice;
		direction[2][2]=-zSlice;
	}
	
	output->SetDirection(direction);
	//---------------------------------------------------------------------	
		
	//Number of images (DWIs+baselines) in the volume
	int nSliceInVolume;
	// figure out how many slices are there in a volume, each unique
	// SliceLocation represent one slice 
	nSliceInVolume = sliceLocations.size();
	int nVolume=nSlice/nSliceInVolume;
	output->SetNumImages(nVolume);
	//------------------------------------------------------------------
	
					
	//Copy MetaDataDictionary from instantiated reader to output image.
	output->SetMetaDataDictionary(gdcmIO->GetMetaDataDictionary());
	this->SetMetaDataDictionary(gdcmIO->GetMetaDataDictionary());
	
	// Allocate the VectorImage
	typedef typename T::IndexType   IndexType;
	IndexType start;
	start.Fill(0);

	VectorImageRegionType region;
	region.SetSize(dimSize);
	region.SetIndex(start);
	
	output->SetLargestPossibleRegion(region);
	output->SetVectorLength(nVolume);
	output->Allocate();
	

  }








template <class T >
void DICOMtoDWIReader<T>::GenerateData( ) {

	typename T::Pointer output = this->GetOutput();

	// allocate the output buffer
	output->SetBufferedRegion( output->GetRequestedRegion() );
	output->Allocate();
		
	InputNamesGeneratorType::Pointer inputNames = InputNamesGeneratorType::New();
	inputNames->SetInputDirectory( m_DirectoryName ); 
 
 	typedef itk::GDCMImageIO ImageIOType;
	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	
	ImageReaderType::Pointer reader = ImageReaderType::New();
	reader->SetImageIO( gdcmIO );

	reader->SetFileNames( inputNames->GetInputFileNames() );
	try
	{
      reader->Update();
    }
	catch (itk::ExceptionObject &excp)
    {
      std::cerr << "Exception thrown while reading the series" << std::endl;
      std::cerr << excp << std::endl;
      return;
    }

	ImageReaderType::DictionaryArrayRawPointer inputDict = reader->GetMetaDataDictionaryArray();
	std::string tag;

	SizeType size=output->GetLargestPossibleRegion().GetSize();
    int nSlice= (output->GetNumImages())*size[2];
		
	// Obtain gradient directions:
	DiffusionDirectionsType DiffusionVectors;
	itk::Array<double> bValues( output->GetNumImages() );
	unsigned int pos = 0;
	tag.clear();
	itk::ExposeMetaData<std::string> ( *(*inputDict)[0], "0018|1312",  tag);
	std::cout<<"row or col"<<tag.c_str()<<std::endl;

	bool xyFlip=false;
	if(strcmp(tag.c_str(), "ROW ")==0){
		xyFlip=true;
	}
	
	for (int k = 0; k < nSlice; k += size[2])
	{
		tag.clear();
		bool exist = itk::ExposeMetaData<std::string> ( *(*inputDict)[k], "0043|1039",  tag);
		float b = atof( tag.c_str() );
		bValues[pos++] = b;
        
		itk::Vector<float,3> vect3d;
		if (!exist || b == 0)
		{
			vect3d.Fill( 0 );
			DiffusionVectors.push_back(vect3d);        
			continue;
		}

		vect3d.Fill( 0 );
		tag.clear();
		itk::ExposeMetaData<std::string> ( *(*inputDict)[k], "0019|10bb",  tag);
		if(xyFlip){
			vect3d[1] = atof( tag.c_str() );
		}else{
			vect3d[0] = atof( tag.c_str() );
		}
		
		tag.clear();
		itk::ExposeMetaData<std::string> ( *(*inputDict)[k], "0019|10bc",  tag);
		if(xyFlip){
			vect3d[0] = atof( tag.c_str() );
		}else{
			vect3d[1] = atof( tag.c_str() );
		}
		
		tag.clear();
		itk::ExposeMetaData<std::string> ( *(*inputDict)[k], "0019|10bd",  tag);
		vect3d[2] = atof( tag.c_str() );

		vect3d.Normalize();
		DiffusionVectors.push_back(vect3d);      
	}

	output->SetBValues(bValues);
  
  // transform gradient directions into RAS frame 
	for (unsigned int k = 0; k < output->GetNumImages(); k++)
    {
          DiffusionVectors[k][2] = -DiffusionVectors[k][2];  // I -> S
    }

	output->CreateDiffusionDirections(DiffusionVectors);
	
	ImageType::Pointer auxImage=ImageType::New();
	auxImage->SetRequestedRegion(output->GetLargestPossibleRegion());
	auxImage=reader->GetOutput();

	typedef itk::ImageRegionIterator< ImageType> ImagesIteratorType;
	ImagesIteratorType it2(auxImage, auxImage->GetRequestedRegion());

	VariableLengthVector<float> aux;
	aux.SetSize(output->GetNumImages());
	
	IndexType index, v_index;

	unsigned int n;	
	
	for ( it2.GoToBegin(); !it2.IsAtEnd(); ++it2 ){
			index=it2.GetIndex();
		
			v_index=index;
			n = index[2]/size[2];
		
			v_index[2]=index[2]-n*size[2];
		
			if(!SliceOrderIS)
			{
				v_index[2]=size[2]-v_index[2]-1;
			}
			
			aux=output->GetPixel(v_index);
			aux[n]=it2.Get();
			output->SetPixel(v_index, aux);
	}
	
}  


} // end namespace itk

#endif
