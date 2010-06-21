/*=========================================================================
 
        Program:        Visualization Toolkit
        Module:         $RCSfile: vtkVolusonReader.h,v $
        Language:       C++
        Date:           $Date: 2001/06/28 10:45:14 $
        Version:        $Version 1.0 $
 
        Autor:          Raul San Jose
                        ETSI Telecomunicacion, Universidad de Valladolid
                        Campus Miguel Delibes, Camino del Cementerio, s/n
                        e-mail: rjosest@atenea.tel.uva.es
 
=========================================================================*/
 
// .NAME vtkVolusonReader - read voluson 530/730 cartesian data file
// .SECTION Description
// vtkVolusonReader is a source object that reads binary voluson 530/730 data files in
// cartesian coordinates and encapsulates the volume in a vtkStructuredPoints data type.

#ifndef __vtkVolusonReader_h
#define __vtkVolusonReader_h

#include "vtkStructuredPointsSource.h"

class VTK_EXPORT vtkVolusonReader : public vtkStructuredPointsSource
{
public:
	vtkVolusonReader();
	~vtkVolusonReader();
	static vtkVolusonReader *New()
		{return new vtkVolusonReader;}
	const char *GetClassName() {return "vtkVolusonReader";};
	void PrintSelf(ostream& os, vtkIndent indent); 

	// Description:
  // Return the MTime also considering the vtkDataReader ivar.
  unsigned long GetMTime();
 
  // Description:
  // Specify file name of voluson data file to read.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  //void SetFileName(char *name);
  //char *GetFileName();
  
  protected:
  void Execute();
  char *FileName;

};
 
#endif
