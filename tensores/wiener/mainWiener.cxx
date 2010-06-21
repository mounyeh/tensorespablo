#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkWienerFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


int main( int argc, char *argv[] ){
  int PNG = 1;
  if( argc!= 5 ){
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "WienerFilter InputImageFile OutputFile iterations lambda" << std::endl;
    return 1;
  } 

  typedef unsigned char   InternalPixelType;
  const   unsigned int    Dimension = 2;
  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
  typedef float OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::WienerFilter< InternalImageType, OutputImageType >  WienerFilterType;
  typedef itk::ImageFileReader< InternalImageType > ReaderType;
  typedef itk::ImageFileWriter< InternalImageType  > WriterType;
  typedef itk::RescaleIntensityImageFilter< OutputImageType,
                        InternalImageType >   CastFilterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  CastFilterType::Pointer caster = CastFilterType::New();

  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  int iterations = atoi(argv[3]);
  float lambda = atof(argv[4]);

  WienerFilterType::Pointer WienerFilter = WienerFilterType::New();
  
  WienerFilter->SetInput(reader->GetOutput());
  WienerFilter->Setiterations(iterations);
  WienerFilter->Setlambda(lambda);
  
  //ThresholdFilter->ThresholdBelow(value);

  if (PNG) {
    caster->SetInput(WienerFilter->GetOutput() );  
    writer->SetInput(caster->GetOutput() );
    caster->SetOutputMinimum( 0 );
    caster->SetOutputMaximum( 255 );
  } else {
    //writer->SetInput(ThresholdFilter->GetOutput() );
  }
  
  writer->Update();

}
