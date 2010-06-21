/*=========================================================================

=========================================================================*/
#ifdef _MSC_VER
#pragma warning( disable : 4611 )
#endif 

#include "itkVOLImageIO.h"
#include <stdio.h>
#include <sys/stat.h>


namespace itk
{


VOLImageIO::VOLImageIO()
{
  m_ChFileName = (char*)NULL;
  m_PstVolFile = (stVOLUMEFILE*)NULL;
}

VOLImageIO::~VOLImageIO()
{
  if( m_PstVolFile != NULL )
    vCloseVolumeFile( m_PstVolFile );
}

void VOLImageIO::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

bool VOLImageIO::CanReadFile( const char* file ) 
{   
    if( m_PstVolFile!=NULL ){
        vCloseVolumeFile( m_PstVolFile );
        m_PstVolFile = NULL;
    }
    m_ChFileName = const_cast<char*>( file );
    if (   ( m_PstVolFile = pstOpenVolumeFile(file) ) == NULL   ){
        m_ChFileName = NULL;
        return false;
    }
    m_EVolType = eVolumeType( m_PstVolFile );
    if ( !(m_EVolType & eCARTESIAN) ){
        m_ChFileName = NULL;
        vCloseVolumeFile( m_PstVolFile );
        m_PstVolFile = NULL;
        return false;
    }
    int iNoOfVolumes = 0;
    if ( m_EVolType & eCFM )
        iNoOfVolumes = 3;
    else if ( m_EVolType & eANGIO )
        iNoOfVolumes = 2;
    else 
        iNoOfVolumes = 1;
    if( iNoOfVolumes<1 ){
        m_ChFileName = NULL;
        vCloseVolumeFile( m_PstVolFile );
        m_PstVolFile = NULL;
        return false;
    }
    if( iNoOfVolumes>1 )
        std::cerr << "Hay mas de un volumen en el fichero; se devolvera solo el primero de ellos" << std::endl;
    return true;
}
  
void VOLImageIO::Read( void* buffer )
{
    if( m_ChFileName==NULL )
       m_ChFileName = const_cast<char*>( m_FileName.c_str() );
    if( m_PstVolFile==NULL ){
       if( !this->CanReadFile( m_ChFileName ) )
           itkExceptionMacro( << "No se pudo leer el fichero" );
    }
    if( m_PstVolFile->pstVolume[0].pbyData==NULL ){
        if( !iLoadVolumeFile(m_PstVolFile) ){
            vCloseVolumeFile( m_PstVolFile );
            itkExceptionMacro( << "No se pudo cargar el volumen" );
        }
    }
    
    memcpy( (BYTE*)buffer, m_PstVolFile->pstVolume[0].pbyData, m_PstVolFile->pstVolume[0].dwDataLength );
}

void VOLImageIO::ReadImageInformation()
{
    if( m_ChFileName==NULL )
       m_ChFileName = const_cast<char*>( m_FileName.c_str() );
    if( m_PstVolFile==NULL ){
       if( !this->CanReadFile( m_ChFileName ) )
           itkExceptionMacro( << "No se pudo leer el fichero" );
    }
    if( m_PstVolFile->pstVolume[0].pbyData==NULL ){
        if( !iLoadVolumeFile(m_PstVolFile) ){
            vCloseVolumeFile( m_PstVolFile );
            itkExceptionMacro( << "No se pudo cargar el volumen" );
        }
    }
    this->SetPixelType( SCALAR );
    this->SetComponentType( UCHAR );
    this->SetNumberOfDimensions( 3 );
    for(int i=0; i<3; i++){
       this->SetSpacing(    i, m_PstVolFile->pstVolume[0].dVoxelSizeMeters );
       this->SetOrigin(     i, 0 );
    }
    this->SetDimensions( 0, m_PstVolFile->pstVolume[0].wSamples );
    this->SetDimensions( 1, m_PstVolFile->pstVolume[0].wShots   );
    this->SetDimensions( 2, m_PstVolFile->pstVolume[0].wImages  );
}  

bool VOLImageIO::CanWriteFile( const char * name )
{
    return false;
}


void VOLImageIO::WriteImageInformation(void)
{
}

void VOLImageIO::Write(const void* buffer)
{
}



} // end namespace itk

