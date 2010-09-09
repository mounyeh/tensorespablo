/*=========================================================================

  Program:   Viewer3D.h
  Language:  C++
  Date:      7-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#ifndef _Viewer3D_h
#define _Viewer3D_h

#include "Fl_VTK_Window.H"
#include <vtkPolyDataMapper.h>
#include <vtkImagePlaneWidget.h>
#include <vtkStreamLine.h>

class Viewer3D : public Fl_VTK_Window
{

public:
  Viewer3D(int x, int y, int w, int h, const char *l);
  virtual ~Viewer3D();
  //void ViewModel(vtkPolyData* data);
  void ConnectMapper(vtkPolyData* data, vtkActor* actor);
  void ConnectMapperFiber(vtkPolyData* data, vtkActor* tractsActor);
  void ShowModel(vtkActor* actor);
  void HideModel(vtkActor* actor);
  void DeleteModel(vtkPolyData* data,vtkActor* actor);
  void Update();
  void SelectSlice(int orientation, unsigned int value);
  void HideSlice(int orientation);
  void SetSliceOpacity(float opacity);
  void SetSliceColor();
  void renderSlice(vtkImageData *image, int orientation);
  void renderSlice(vtkImageData *image, int orientation, int slice);
  void renderTractsTube(vtkPolyData *tracts, vtkActor* tractsActor);
  void renderTracts(vtkPolyData *tracts, vtkActor* tractsActor, double color[3], double radius);
  void renderTractsNew(vtkStreamLine *tracts, vtkActor* tractsActor);
  void ChangeColor(vtkActor* actor, double r, double g, double b);
  void SetScalarRange(vtkActor* actor,double min,double max);
  void ChangeOpacity(vtkActor* actor, float value);
  void ChangeSpecular(vtkActor* actor, float value);
  void ChangeDiffuse(vtkActor* actor, float value);
  void ChangeAmbient(vtkActor* actor, float value);
  void ChangeSpecularPower(vtkActor* actor, float value);
  void SetSliceColorDefault();
  void ResetView(void);
  void SetZoom(float value);
  void SetView(unsigned int value);
  void Load();
  
  vtkImagePlaneWidget* planeWidgetX;
  vtkImagePlaneWidget* planeWidgetY;
  vtkImagePlaneWidget* planeWidgetZ;
private:
  const char * m_filename;
  
};
#endif




