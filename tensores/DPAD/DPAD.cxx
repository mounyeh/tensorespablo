#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDPADFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkComputeDiffusionTermFilter.h"
#include "itkCommand.h"


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
		std::cerr << argv[0] << " input output [ iterations=5 tau=2 radius=2 useMode=0 useAOS=1 diffusionTermMode=[Aja,Yu,Simplified] ]" << std::endl;
		return 1;
	}
	unsigned int iterations  = 5;
	float        tau         = 2.0f;
	unsigned int radius      = 2;
	bool         mode        = false;
	bool         aos         = true;
	if( argc>3 )
		iterations  = (unsigned int)( ::atoi(argv[3]) );
	if( argc>4 )
		tau         = (float)( ::atof(argv[4]) );
	if( argc>5 )
		radius      = (unsigned int)( ::atoi(argv[5]) );
	if( argc>6 )
		mode        = (bool)( ::atoi(argv[6]) );
	if( argc>7 )
		aos         = (bool)( ::atoi(argv[7]) );
	int dtM = 0;
	if( argc>8 ){
		if( !strcmp(argv[8],"Aja") )
			dtM = 0;
		else if( !strcmp(argv[8],"Yu") )
			dtM = 1;
		else if( !strcmp(argv[8],"Simplified") )
			dtM = 2;
		else{
			std::cerr << "??? Unknown way to compute the diffusion term from the coefficient of variation" << std::endl;
			exit( -1 );
		}
	}
	
	//--------------------------------------------------------------------------------------------------------------------------
	
	typedef   float          InputPixelType;
	typedef   unsigned char  OutputPixelType;
	
	typedef itk::Image< InputPixelType,  2 >   InputImageType;
	typedef itk::Image< OutputPixelType, 2 >   OutputImageType;
	
	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	
	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();
	
	reader->SetFileName( argv[1] );
	writer->SetFileName( argv[2] );
	
	typedef itk::DPADFilter< InputImageType, InputImageType >                    FilterType;
	typedef itk::RescaleIntensityImageFilter<  InputImageType, OutputImageType > RescaleFilterType;
	typedef FilterType::DiffusionTermType                                        DiffusionType;

	FilterType::Pointer         filter        = FilterType::New();
	DiffusionType::Pointer      diffusion     = DiffusionType::New();
	RescaleFilterType::Pointer  rescaleFilter = RescaleFilterType::New();
	FilterType::InputSizeType   rad;
	rad.Fill( radius );
	
	filter->SetInput( reader->GetOutput() );
	filter->SetIter( iterations );
	filter->SetTimeStep( tau );
	filter->SetRadius( rad );
	filter->SetUseMode( mode );
	filter->SetUseAOS( aos );
	switch( dtM ){
	case 0:
		filter->SetUseAjaDiffusionTerm();
		break;
	case 1:
		filter->SetUseYuDiffusionTerm();
		break;
	case 2:
		filter->SetUseSimplifiedDiffusionTerm();
		break;
	}
	filter->SetDiffusionTerm( diffusion );
	
	
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

