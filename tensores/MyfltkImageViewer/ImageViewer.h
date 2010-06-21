/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageViewer.h,v $
  Language:  C++
  Date:      $Date: 2003/07/03 16:03:08 $
  Version:   $Revision: 1.2 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __ImageViewer_h
#define __ImageViewer_h

#include "ImageViewerGUI.h"

//namespace fltk {

template <class ImagePixelType, class OverlayPixelType>
class ImageViewer : public ImageViewerGUI
{
public:

   /**
   * Standard "Self" typedef.
   */
  typedef ImageViewer         Self;

  /** 
   * Smart pointer typedef support.
   */
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

 
  typedef itk::Image< ImagePixelType, 3 >   ImageType;
  typedef itk::Image< OverlayPixelType, 3 > OverlayType;
  typedef GLSliceView< ImagePixelType, OverlayPixelType > GLSliceViewType;
  typedef typename GLSliceViewType::ColorTablePointer ColorTablePointer;
  
  ImageViewer();
  virtual ~ImageViewer();
  virtual void SetImage(ImageBase<3> * img);
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
  virtual void SetOrientation(void);
  virtual void CreateGLSliceView( Fl_Group * g , Fl_Gl_Window * w );
  virtual void SetLabel( const char * );
  virtual void ShowClickedPointsWindow(void);
  virtual void UpdateClickedPoints(void);
  virtual void ClearClickedPoints(void);
  virtual void ClickSelectCallBack(
                void (*newClickSelectArgCallBack)(float, float,
                                                  float, float,int,
                                                  void *),
                     void * newClickSelectArg);
  
  virtual GLSliceView<ImagePixelType, OverlayPixelType>* GetSliceViewPointer();
 
  virtual bool viewDetails(void);
  virtual bool viewValue(void);
  virtual bool viewCrosshairs(void);
  virtual bool flipX(void);
  virtual bool flipY(void);
  virtual bool flipZ(void);

  virtual void ViewDetails(bool detail);
  virtual void ViewValue(bool value);
  virtual void ViewCrosshairs(bool crosshairs);
  virtual void flipX(bool t);
  virtual void flipY(bool t);
  virtual void flipZ(bool t);

  /** Specify the opacity of the overlay */
  virtual void  SetOverlayOpacity(float newOverlayOpacity);
  
  /** Get the opacity of the overlay */
  virtual float GetOverlayOpacity(void) const;
  
  /** Show slider to control overlay opacity */
  virtual void ShowOverlayOpacityControl(void);
  
  /** Get the ColorTable for the Overlay */
  virtual ColorTablePointer GetOverlayColorTable(void);

  /** Set the ColorTable for the Overlay */
  virtual void SetOverlayColorTable(ColorTablePointer newColorTable);


private:
  GLSliceViewType * glSliceView;
};


//} // end namespace fltk

#ifndef ITK_MANUAL_INSTANTIATION
#include <ImageViewer.txx>
#endif




#endif
