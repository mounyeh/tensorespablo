#include "itkMLXB.h"

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
		std::cerr << "Registration level:  " << m_Level << std::endl;
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
		std::cerr << "Usage: " << argv[0] << "   FixedImagefile   MovingImageFile   ";
		std::cerr << "RegisteredImageFile [DownSample=6] [BlockSize=2] [SearchSize=1] [NLevels=5] [GWidth=4] [Sigma=16.0] [NIter=5] [Metric=6] [Bins=20] [Lambda=0.5]" << std::endl;
		return -1;
	}

	const    unsigned int    Dimension = 3;
	typedef  unsigned char   PixelType;

	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;
	typedef MovingImageType::PointType          PointType;

	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	typedef itk::MLXB< FixedImageType, MovingImageType > RegisterType;

	RegisterType::Pointer reg = RegisterType::New();

	typedef AlgorithmInterfaceCommand<RegisterType>     CommandType;
	CommandType::Pointer command = CommandType::New();
    reg->AddObserver( itk::IterationEvent(), command );

	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	
	fixedImageReader->SetFileName(  argv[1] );
	movingImageReader->SetFileName( argv[2] );
	
	fixedImageReader->Update( );
	movingImageReader->Update();

	if( argc >4 ){
		FixedImageType::SizeType sample;
		sample.Fill( ::atoi(argv[4]) );
		reg->SetSample( sample );
	}
	if( argc >5 ){
		FixedImageType::SizeType block;
		block.Fill( ::atoi(argv[5]) );
		reg->SetBlockSize( block );
	}
	if( argc >6 ){
		FixedImageType::SizeType search;
		search.Fill( ::atoi(argv[6]) );
		reg->SetSearchSize( search );
	}
	if( argc >7 ){
		reg->SetNLevels( ::atoi(argv[7]) );
	}
	if( argc >8 ){
		RegisterType::SizeType size;
		size.Fill( ::atoi(argv[8]) );
		reg->SetGaussianRadius( size );
	}

	itk::Array<double> sigma;
	double _sigma = 16.0;
	if( argc>9 )
		_sigma = (double)( ::atof(argv[9]) );
	sigma.SetSize(3);
	sigma[0] = _sigma;
	sigma[1] = _sigma;
	sigma[2] = _sigma;

	reg->SetSigma( sigma );

	if( argc>10 )
		reg->SetNIter( ::atoi(argv[10]) );

	if( argc>11 )
		reg->SetMetric( ::atoi(argv[11]) );
	switch( reg->GetMetric() ){
		case 1:
			std::cerr << "Using NMI" << std::endl;
			break;
		case 2:
			std::cerr << "Using MI" << std::endl;
			break;
		default:
			std::cerr << "Using pondered NMI and MI" << std::endl;
			break;
	}

	if( argc>12 )
		reg->SetBins(   (unsigned int)( ::atoi(argv[12]) )   );

	if( argc>13 )
		reg->SetLambda(   (double)( ::atof(argv[13]) )   );

	reg->SetInput1( fixedImageReader->GetOutput() );
	reg->SetInput2(	movingImageReader->GetOutput() );
	reg->Start();

	typedef itk::ImageFileWriter< MovingImageType >  WriterType;
	WriterType::Pointer      writer =  WriterType::New();

	writer->SetFileName( argv[3] );
	writer->SetInput( reg->GetOutput()   );
	
	try{
		writer->Update();
	}
	catch( itk::ExceptionObject & err ) 
    {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return -1;
    }
}
