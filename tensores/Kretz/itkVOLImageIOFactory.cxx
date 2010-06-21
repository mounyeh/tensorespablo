/*=========================================================================
=========================================================================*/
#include "itkVOLImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkVOLImageIO.h"
#include "itkVersion.h"

namespace itk
{

VOLImageIOFactory::VOLImageIOFactory()
{
  this->RegisterOverride("itkImageIOBase",
                         "itkVOLImageIO",
                         "VOL Image IO",
                         1,
                         CreateObjectFunction<VOLImageIO>::New());
}
  
VOLImageIOFactory::~VOLImageIOFactory()
{
}

const char* 
VOLImageIOFactory::GetITKSourceVersion(void) const
{
  return ITK_SOURCE_VERSION;
}

const char* 
VOLImageIOFactory::GetDescription(void) const
{
  return "This allows to read GE Kretz files in cartesian coordinates";
}

} // end namespace itk

