#include "itkImage.h"
#include "itkImageFileReader.h"

int main( int argc, char * argv[] )
{
  if( argc < 2 ){
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile1 inputImageFile2 ... inputImageFileN" << std::endl;
    return 1;
  }
  
  typedef   unsigned char  InputPixelType;
  typedef itk::Image< InputPixelType,  3 >   InputImageType;
  typedef itk::ImageFileReader< InputImageType  >   ReaderType;
  
  ReaderType::Pointer reader = ReaderType::New();
  unsigned int _MAX = 0;
  for( unsigned int k=1; k<argc;k++ ){
      if( strlen(argv[k])>_MAX )
          _MAX = ::strlen(argv[k]);
  }
  _MAX++;
  for( unsigned int k=1; k<argc; k++ ){
      reader->SetFileName( argv[k] );
      unsigned int ll    = ::strlen( argv[k] );
      unsigned int blank = _MAX - ll;
      try{
          reader->Update();
          std::cout << argv[k];
          for( unsigned int l=0; l<blank; l++ )
              std::cout << " ";
          std::cout << ":  Se puede leer" << std::endl;
      }
      catch( itk::ExceptionObject & err ){
          std::cout << argv[k];
          for( unsigned int l=0; l<blank; l++ )
              std::cout << " ";
          std::cout << ":  No se puede leer !!!" << std::endl;
      }
  }
  
  
  return 0;
}
