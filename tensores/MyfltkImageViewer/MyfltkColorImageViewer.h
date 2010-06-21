/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: fltkImageViewer.h,v $
  Language:  C++
  Date:      $Date: 2005/06/13 13:45:08 $
  Version:   $Revision: 1.12 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __MyfltkColorImageViewer_h
#define __MyfltkColorImageViewer_h

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Choice.H>
#include <GLColorSliceView.h>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Browser.H>

namespace fltk {

template <class ImagePixelType, class OverlayPixelType>
  class MyColorImageViewer : public GLColorSliceView<ImagePixelType, OverlayPixelType>
{
public:

   /**
   * Standard "Self" typedef.
   */
  typedef MyColorImageViewer         Self;

  /** 
   * Smart pointer typedef support.
   */
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

 
  typedef itk::Image< ImagePixelType, 3 >   ImageType;
  typedef itk::Image< OverlayPixelType, 3 > OverlayType;
  typedef GLColorSliceView< ImagePixelType, OverlayPixelType > GLColorSliceViewType;
  typedef typename GLColorSliceViewType::ColorTablePointer ColorTablePointer;
    
  MyColorImageViewer(int x,int y,int w,int h, const char * label=0);
  virtual ~MyColorImageViewer();
  //virtual void SetImage(ImageBase<3> * img);
  virtual void SetOverlay(ImageBase<3> * img);
  virtual void Show(void);
  virtual void Hide(void);
  virtual void Update(void);
  virtual void Synchronize(void);
  virtual void ImageMode(ImageModeType mode);
  virtual void SelectSlice(unsigned int num);
  virtual void SetIntensityWindowingMin(float val);
  virtual void SetIntensityWindowingMax(float val);
  virtual void CenterWindow(void);
  virtual void ZoomIn(void);
  virtual void ZoomOut(void);
  virtual void ShiftUp(void);
  virtual void ShiftDown(void);
  virtual void ShiftLeft(void);
  virtual void ShiftRight(void);
  virtual void SetOrientation(int val);
  virtual void ClickSelectCallBack(
                void (*newClickSelectArgCallBack)(float, float,
                                                  float, float,int,
                                                  void *),
                     void * newClickSelectArg);

  virtual void ChangeActiveClass(int value);   
  virtual void ClearClickedPoints(void);
  virtual void ViewDetails(bool detail);
  virtual void ViewValue(bool value);
  virtual void ViewCrosshairs(bool crosshairs);
  virtual unsigned int GetNumSlices();
  /** Specify the opacity of the overlay */
  virtual void  SetOverlayOpacity(float newOverlayOpacity);
  virtual int  handle(int event);


private:
  //GLSliceViewType * glSliceView;
  


};


} // end namespace fltk

#ifndef ITK_MANUAL_INSTANTIATION
#include "MyfltkColorImageViewer.txx"
#endif




#endif
