/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkNrrdVectorImageReadWriteTest.cxx,v $
 Language:  C++
 Date:      $Date: 2005/08/20 22:47:07 $
 Version:   $Revision: 1.1 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
/*
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
*/
#include <fstream>
#include "itkDICOMtoDWIReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
//#include "itkVectorImage.h"
#include "itkDWImages.h"
#include "itkLMMSEVectorImageFilter.h"
#include "itkDTIEstimateTensorFilter.h"
#include "itkVector.h"
#include <string.h>
#include "itkDTITensor.h"
#include "itkComputeFAFilter.h"

int main( int ac, char* av[] )
{
	if(ac < 3)
    {
		std::cerr << "Usage: " << av[0] << " Input Output [radiusLMMSEest=1 radiusLMMSEfilt=2 iterLMMSE=1 iterEstim=1 ]\n";
		return EXIT_FAILURE;
    }

	typedef  float                                           PixelType;
	typedef itk::DWImages<PixelType, 3>                      ImageType;
	typedef itk::LMMSEVectorImageFilter<ImageType,ImageType> FilterType;
	
	itk::DICOMtoDWIReader<ImageType>::Pointer reader = itk::DICOMtoDWIReader<ImageType>::New();
	reader->SetDirectoryName(av[1]);
	
	
	try
    {
		reader->Update();
    }
	catch ( itk::ExceptionObject & e )
    {
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e.GetDescription() << std::endl;
		std::cerr << e.GetLocation() << std::endl;
		return EXIT_FAILURE;
    }
	
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( reader->GetOutput() );
	
	//==========================================================================
	// Set the parameters to the filter:
	FilterType::InputSizeType radius;
	unsigned int rade = 1;
	if( ac>3 )
		rade = ::atoi( av[3] );
	radius.Fill( rade );
//-----------------------
	//radius[2] = 0;
//-----------------------
	filter->SetRadiusEstimation( radius );
	unsigned int radf = 2;
	if( ac>4 )
                radf = ::atoi( av[4] );
	radius.Fill( radf );
	radius[2] = 0;
	filter->SetRadiusFiltering( radius );
	unsigned int iter1 = 1;
	if( ac>5 )
		iter1 = ::atoi( av[5] );
	
	filter->SetIterations( iter1 );
	filter->SetUseAbsoluteValue( false );
	filter->SetKeepValue( true );
	filter->SetMinimumNumberOfUsedVoxelsEstimation( 5 );
	filter->SetMinimumNumberOfUsedVoxelsFiltering( 5 );
	filter->SetMinimumNoiseSTD( 0 );
	filter->SetMaximumNoiseSTD( 32000 );
	filter->SetFirstBaseline( 0 );
	filter->SetHistogramResolutionFactor( 2.0 );
	filter->SetChannels( reader->GetOutput()->GetVectorLength() );
	//==========================================================================
	
	try
    {
		if( iter1>0 ){
			std::cerr << "Filtrando..." << std::endl;
			filter->Update();
			std::cerr << "Filtrado!" << std::endl;
	std::cerr << "===========================================================" << std::endl;
        std::cerr << "                  CONFIGURE DONE!!!!"                        << std::endl;
        std::cerr << "Channels:    " << reader->GetOutput()->GetNumDWImages()         << std::endl;
        std::cerr << "NBaselines:  " << reader->GetOutput()->GetNumBaselines()        << std::endl;
        std::cerr << "DWIChannels: " << reader->GetOutput()->GetIndexesDWImages()     << std::endl;
        std::cerr << "Baselines:   " << reader->GetOutput()->GetIndexesBaselines()    << std::endl;                                                                                                                                                                              
        std::cerr << "B:           " << reader->GetOutput()->GetB_Value()              << std::endl;
	std::cerr << "===========================================================" << std::endl;
			std::cerr << "Copiando paramteros..." << std::endl;
			filter->GetOutput()->SetNumImages( reader->GetOutput()->GetNumImages() );
			filter->GetOutput()->SetNumDWImages( reader->GetOutput()->GetNumDWImages() );
			filter->GetOutput()->SetNumBaselines( reader->GetOutput()->GetNumBaselines() );
			filter->GetOutput()->SetIndexesDWImages( reader->GetOutput()->GetIndexesDWImages() );
			filter->GetOutput()->SetIndexesBaselines( reader->GetOutput()->GetIndexesBaselines() );
			filter->GetOutput()->SetBValues( reader->GetOutput()->GetBValues() );
			ImageType::DiffusionDirectionsType directions;
			reader->GetOutput()->GetDiffusionDirections( directions );
			filter->GetOutput()->SetDiffusionDirections( directions );
	std::cerr << "===========================================================" << std::endl;                                                                                                                                                                 
        std::cerr << "                  CONFIGURE DONE!!!!"                        << std::endl;                                                                                                                                                                 
        std::cerr << "Channels:    " << filter->GetOutput()->GetNumDWImages()         << std::endl;                                                                                                                                                              
        std::cerr << "NBaselines:  " << filter->GetOutput()->GetNumBaselines()        << std::endl;                                                                                                                                                              
        std::cerr << "DWIChannels: " << filter->GetOutput()->GetIndexesDWImages()     << std::endl;                                                                                                                                                              
        std::cerr << "Baselines:   " << filter->GetOutput()->GetIndexesBaselines()    << std::endl;                                                                                                                                                                              
        std::cerr << "B:           " << filter->GetOutput()->GetB_Value()              << std::endl;                                                                                                                                                              
        std::cerr << "===========================================================" << std::endl;
			std::cerr << "Copiados" << std::endl;
		}
    }
	catch ( itk::ExceptionObject & e )
    {
		std::cerr << "exception inL MMSE filter" << std::endl;
		std::cerr << e.GetDescription() << std::endl;
		std::cerr << e.GetLocation() << std::endl;
		return EXIT_FAILURE;
    }



	typedef itk::DTITensor<float>                              TPixelType;	
	typedef itk::Image<TPixelType,3>                           TImageType;
	typedef itk::DTIEstimateTensorFilter<ImageType,TImageType> DTIType;
	DTIType::Pointer dti = DTIType::New();
	if( iter1>0 )
		dti->SetInput( filter->GetOutput() );
	else
		dti->SetInput( reader->GetOutput() );
	// The number of iterations; by default, only one iteration (simple WLS) is
	// performed:
	unsigned int iter2 = 1;
	if( ac>6 )
		iter2 = ::atoi( av[6] );
	dti->SetIterations( iter2 );
	std::cerr << "Configuring DTI estimator..." << std::endl;
	dti->Configure();
	std::cerr << "Done." << std::endl;
	
	std::cout << "Computing mask..." << std::endl;
	dti->ComputeMask();
	std::cout << "Done. Threshold: " << dti->GetThreshold() << std::endl;
	//==========================================================================	
	
	
	
	
	try
    {
		std::cout << "Computing the tensor..." << std::endl;
		dti->Update();
		//==================================================================================================
		typedef itk::Image<unsigned char,3>    MaskType;
		typedef itk::ImageFileWriter<MaskType> MWType;
		MWType::Pointer mw = MWType::New();
		mw->SetFileName( "Tensor_Mask_Computed.mhd" );
		mw->SetInput( dti->GetMask() );
		mw->Update();
		//==================================================================================================
		std::cout << "Done" << std::endl;
    }
	catch ( itk::ExceptionObject & e )
    {
		std::cerr << "exception in file writer " << std::endl;
		std::cerr << e.GetDescription() << std::endl;
		std::cerr << e.GetLocation() << std::endl;
		return EXIT_FAILURE;
    }

	typedef float                                       SPixelType;
	typedef itk::Image<SPixelType,3>                    SImageType;
	typedef itk::ComputeFAFilter<TImageType,SImageType> FAType;	
	FAType::Pointer fa = FAType::New();
	fa->SetInput( dti->GetOutput() );
	fa->Update();
	
	// Generate test image
	itk::ImageFileWriter<SImageType>::Pointer writer = itk::ImageFileWriter<SImageType>::New();
	writer->SetInput( fa->GetOutput() );
	writer->SetFileName(av[2]);
	try
    {
		writer->Update();
    }
	catch ( itk::ExceptionObject & e )
    {
		std::cerr << "exception in file writer " << std::endl;
		std::cerr << e.GetDescription() << std::endl;
		std::cerr << e.GetLocation() << std::endl;
		return EXIT_FAILURE;
    }
	return EXIT_SUCCESS;
}
