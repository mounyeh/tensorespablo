/*=========================================================================

=========================================================================*/
#ifndef __itkVOLImageIO_h
#define __itkVOLImageIO_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <fstream>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "itkImageIOBase.h"
#include "importvol.h"

#define VOLUME_NO(wTag) ((wTag & 0x0F00) >> 8)
#define FILE_ID "KRETZFILE 1.0   "

namespace itk
{





class ITK_EXPORT VOLImageIO : public ImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef VOLImageIO          Self;
  typedef ImageIOBase         Superclass;
  typedef SmartPointer<Self>  Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VOLImageIO, ImageIOBase);

  /*-------- This part of the interface deals with reading data. ------ */

  virtual bool SupportsDimension( unsigned long dim )
  {return(dim==3);}

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanReadFile(const char*);
  
  /** Set the spacing and diemention information for the set filename. */
  virtual void ReadImageInformation();
  
  /** Reads the data from disk into the memory buffer provided. */
  virtual void Read(void* buffer);

  /*-------- This part of the interfaces deals with writing data. ----- */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanWriteFile(const char*);

  /** Writes the spacing and dimentions of the image.
   * Assumes SetFileName has been called with a valid file name. */
  virtual void WriteImageInformation();

  /** Writes the data to disk from the memory buffer provided. Make sure
   * that the IORegion has been set properly. */
  virtual void Write(const void* buffer);
protected:
  VOLImageIO();
  ~VOLImageIO();
  void PrintSelf( std::ostream& os, Indent indent ) const;
private:
  VOLImageIO(const Self&);     //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  char*         m_ChFileName;
  stVOLUMEFILE* m_PstVolFile;
  eVOLUMETYPE   m_EVolType;
};

} // end namespace itk

#endif // __itkVOLImageIO_h

