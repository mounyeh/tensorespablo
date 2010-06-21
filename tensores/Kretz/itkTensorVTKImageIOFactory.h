/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorVTKImageIOFactory.h,v $
  Language:  C++
  Date:      $Date: 2007/03/22 14:28:53 $
  Version:   $Revision: 1.10 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorVTKImageIOFactory_h
#define __itkTensorVTKImageIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
/** \class TensorVTKImageIOFactory
 * \brief Create instances of VTKImageIO objects using an object factory.
 */
class ITK_EXPORT TensorVTKImageIOFactory : public ObjectFactoryBase
{
public:  
  /** Standard class typedefs. */
  typedef TensorVTKImageIOFactory        Self;
  typedef ObjectFactoryBase        Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Class Methods used to interface with the registered factories. */
  virtual const char* GetITKSourceVersion(void) const;
  virtual const char* GetDescription(void) const;
    
  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorVTKImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
    {
    TensorVTKImageIOFactory::Pointer vtkFactory = TensorVTKImageIOFactory::New();
    ObjectFactoryBase::RegisterFactory(vtkFactory);
    }
  
protected:
  TensorVTKImageIOFactory();
  ~TensorVTKImageIOFactory();

private:
  TensorVTKImageIOFactory(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

  
} // end namespace itk

#endif
