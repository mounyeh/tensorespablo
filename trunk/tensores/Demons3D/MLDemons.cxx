#include "itkMLD.h"

#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"

template <typename TAlgorithm>
class AlgorithmInterfaceCommand : public itk::Command 
{
public:
	typedef  AlgorithmInterfaceCommand      Self;
	typedef  itk::Command                   Superclass;
	typedef  itk::SmartPointer<Self>        Pointer;
	itkNewMacro( Self );
protected:
	AlgorithmInterfaceCommand()
	{
		m_Level=0;
	};
public:
	typedef   TAlgorithm                              AlgorithmType;
	typedef   AlgorithmType *                         AlgorithmPointer;
	
	void Execute(itk::Object * object, const itk::EventObject & event)
	{
		if( typeid( event ) != typeid( itk::IterationEvent ) )
			return;
		std::cout << "Registration level:  " << m_Level << std::endl;
		m_Level++;
	}
	void Execute(const itk::Object * , const itk::EventObject & )
	{
		return;
	}
private:
	unsigned int m_Level;
};

int main( int argc, char *argv[] )
{
	if( argc<4 ){
		std::cerr << "Mising parameters!" << std::endl;
		std::cerr << "Usage: " << argv[0] << "   FixedImagefile   MovingImageFile   RegisteredImageFile [NLevels=5] [NSteps=10] [UseEReg=1] [UseFreg=0] [SigmaE=2] [SigmaF=9] [SigmaStats=4.1410] [SigmaGradient=0.5770]" << std::endl;
		return -1;
	}

	const    unsigned int    Dimension = 3;
	typedef  unsigned char   PixelType;
	typedef  unsigned char   OPixelType;

	typedef itk::Image< PixelType, Dimension >          FixedImageType;
	typedef itk::Image< PixelType, Dimension >          MovingImageType;
	typedef itk::Image< OPixelType, Dimension >         OutputImageType;
	typedef MovingImageType::PointType                  PointType;

	typedef itk::ImageFileReader< FixedImageType  >     FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType >     MovingImageReaderType;

	typedef itk::MLD< FixedImageType, MovingImageType > RegisterType;

	RegisterType::Pointer reg = RegisterType::New();

	FixedImageReaderType::Pointer  fixedImageReader  =  FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader =  MovingImageReaderType::New();
	
	fixedImageReader->SetFileName(  argv[1] );
	movingImageReader->SetFileName( argv[2] );
	
	fixedImageReader->Update( );
	movingImageReader->Update();


	typedef AlgorithmInterfaceCommand<RegisterType>     CommandType;
	CommandType::Pointer command = CommandType::New();
    	reg->AddObserver( itk::IterationEvent(), command );

	reg->SetInput1( fixedImageReader->GetOutput() );
	reg->SetInput2(	movingImageReader->GetOutput() );
	if( argc>4 )
		reg->SetNLevels(   ::atoi(argv[4])   );
	if( argc>5 )
		reg->SetSteps(     ::atoi(argv[5])   );
	if( argc>6 )
		reg->SetUseElasticRegularization(   (bool)::atoi(argv[6])   );
	if( argc>7 )
		reg->SetUseFluidRegularization(   (bool)::atoi(argv[7])   );
	if( argc>8 )
		reg->SetSigmaElastic(   (int)::atoi(argv[8])   );
	if( argc>9 )
		reg->SetSigmaFluid(   (int)::atoi(argv[9])   );
	if( argc>10 )
		reg->SetSigmaStats(   (float)::atof(argv[10])   );
	if( argc>11 )
		reg->SetSigmaGradient(   (float)::atof(argv[11])   );
	reg->Start();

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	WriterType::Pointer      writer =  WriterType::New();

	writer->SetFileName( argv[3] );
	writer->SetInput( reg->GetOutput()   );
	
	try{
		writer->Update();
	}
	catch( itk::ExceptionObject & err ) 
    	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return -1;
    	}
}
