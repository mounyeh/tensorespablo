/*=========================================================================
=========================================================================*/
#ifndef __itkVOLImageIOFactory_h
#define __itkVOLImageIOFactory_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{

class ITK_EXPORT VOLImageIOFactory : public ObjectFactoryBase
{
public:  
  /** Standard class typedefs. */
  typedef VOLImageIOFactory         Self;
  typedef ObjectFactoryBase         Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Class methods used to interface with the registered factories. */
  virtual const char* GetITKSourceVersion(void) const;
  virtual const char* GetDescription(void) const;
    
  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);
  static VOLImageIOFactory* FactoryNew() { return new VOLImageIOFactory;}
  /** Run-time type information (and related methods). */
  itkTypeMacro(VOLImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
  {
    VOLImageIOFactory::Pointer VOLFactory = VOLImageIOFactory::New();
    ObjectFactoryBase::RegisterFactory(VOLFactory);
  }
  
protected:
  VOLImageIOFactory();
  ~VOLImageIOFactory();

private:
  VOLImageIOFactory(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
  
  
} // end namespace itk

#endif
