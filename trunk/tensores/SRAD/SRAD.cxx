#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSRADFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSRADDiffusionTensorFilter.h"
#include "itkCommand.h"
#include <time.h>


class FilterInterfaceCommand : public itk::Command 
{
public:
	typedef  FilterInterfaceCommand    Self;
	typedef  itk::Command              Superclass;
	typedef  itk::SmartPointer<Self>   Pointer;
	itkNewMacro( Self );
protected:
	FilterInterfaceCommand(){}
public:	
	void Execute(itk::Object * object, const itk::EventObject & event)
	{
		if( typeid( event ) != typeid( itk::IterationEvent ) )
			return;
		std::cerr << " ";
	}	
	void Execute(const itk::Object * , const itk::EventObject & ){ return; }
};



int main( int argc, char * argv[] )
{
	//--------------------------------------------------------------------------------------------------------------------------
	// Check arguments
	if( argc < 3 ){
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " input output [ expConstant=25.0 iterations=50 tau=0.05 beta=0.02f" << std::endl;
		return 1;
	}
	float        expConstant = 25.0f;
	unsigned int iterations  = 50;
	float        tau         = 0.05f;
	float        beta        = 0.02f;
	if( argc>3 )
		expConstant = (unsigned int)( ::atoi(argv[3]) );
	if( argc>4 )
		iterations  = (unsigned int)( ::atoi(argv[4]) );
	if( argc>5 )
		tau         = (float)( ::atof(argv[5]) );
	if( argc>6 )
		beta        = (float)( ::atof(argv[6]) );
	
	//--------------------------------------------------------------------------------------------------------------------------
	
	typedef   unsigned char  InputPixelType;
	typedef   unsigned char  OutputPixelType;
	
	typedef itk::Image< InputPixelType,  2 >   InputImageType;
	typedef itk::Image< OutputPixelType, 2 >   OutputImageType;
	
	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	
	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();
	
	reader->SetFileName( argv[1] );
	writer->SetFileName( argv[2] );
	
	typedef itk::SRADFilter< InputImageType, InputImageType >                     FilterType;
	typedef itk::RescaleIntensityImageFilter<  InputImageType, OutputImageType > RescaleFilterType;
	typedef FilterType::TensorFilterType                                         DiffusionType;

	FilterType::Pointer         filter        = FilterType::New();
	DiffusionType::Pointer      diffusion     = DiffusionType::New();
	RescaleFilterType::Pointer  rescaleFilter = RescaleFilterType::New();
	
	filter->SetInput( reader->GetOutput() );
	filter->SetExpConstant( expConstant );
	filter->SetIter( iterations );
	filter->SetTimeStep( tau );
	filter->SetTensorFilter( diffusion );
	filter->SetBeta(     beta   );
	
	for( unsigned int k=0; k<iterations; k++ )
		std::cerr << ".";
	for( unsigned int k=0; k<iterations; k++ )
		std::cerr << (char)8;
	
	
	FilterInterfaceCommand::Pointer command = FilterInterfaceCommand::New();
	filter->AddObserver( itk::IterationEvent(), command );
	rescaleFilter->SetInput( filter->GetOutput() );
	writer->SetInput( rescaleFilter->GetOutput() );
	rescaleFilter->SetOutputMinimum(  0  );
	rescaleFilter->SetOutputMaximum( 255 );
	try{
		writer->Update();
		std::cerr << (char)7 << '\r';
	}
	catch( itk::ExceptionObject & err ) 
    {
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return -1;
    }
	
	return 0;
}

