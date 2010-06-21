/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorVTKImageIOFactory.cxx,v $
  Language:  C++
  Date:      $Date: 2007/03/22 14:28:53 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkTensorVTKImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkTensorVTKImageIO.h"
#include "itkVersion.h"

namespace itk
{

TensorVTKImageIOFactory::TensorVTKImageIOFactory()
{
  this->RegisterOverride("itkImageIOBase",
                         "itkTensorVTKImageIO",
                         "VTK Image IO",
                         1,
                         CreateObjectFunction<TensorVTKImageIO>::New());
}
  
TensorVTKImageIOFactory::~TensorVTKImageIOFactory()
{
}

const char* 
TensorVTKImageIOFactory::GetITKSourceVersion(void) const
{
  return ITK_SOURCE_VERSION;
}

const char* 
TensorVTKImageIOFactory::GetDescription(void) const
{
  return "VTK ImageIO Factory, allows the loading of VTK images into ITK";
}

} // end namespace itk
