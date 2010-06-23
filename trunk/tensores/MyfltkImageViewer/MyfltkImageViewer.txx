/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: fltkImageViewer.txx,v $
  Language:  C++
  Date:      $Date: 2003/09/28 15:23:29 $
  Version:   $Revision: 1.13 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _MyfltkImageViewer_txx
#define _MyfltkImageViewer_txx

#include "MyfltkImageViewer/MyfltkImageViewer.h"


namespace fltk {
  
template <class ImagePixelType, class OverlayPixelType>
MyImageViewer<ImagePixelType,OverlayPixelType>
::MyImageViewer(int x,int y,int w,int h, const char * label)
:GLSliceView<ImagePixelType,OverlayPixelType>(x,y,w,h,label) 
{

   	
}


  
template <class ImagePixelType, class OverlayPixelType>
MyImageViewer<ImagePixelType,OverlayPixelType>
::~MyImageViewer()
{

}

template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::SetImage(ImageBase<3> * img)
{
  ImageType * image = dynamic_cast<ImageType *>( img );
  this->SetInputImage( image );
  Synchronize();
}


template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::SetOverlay(ImageBase<3> * img)
{
  OverlayType * overlay = dynamic_cast<OverlayType *>( img );
  this->SetInputOverlay( overlay );
  Synchronize();
}

/* Specify the opacity of the overlay */
template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::SetOverlayOpacity(float newOverlayOpacity)
{
  this->OverlayOpacity( newOverlayOpacity );
  this->update();
}


template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::Show(void)
{
  static bool firstTime = true;
  this->show();
  this->update();
  if( firstTime )
  {
    firstTime = false;
    Fl::check();
    this->redraw();
    Fl::check();
  }
}

template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::Hide(void)
{
  this->hide();
  this->update();
  Fl::check();
}


template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::Update(void)
{
  this->update();
}

template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::ImageMode(ImageModeType mode)
{
  this->imageMode(mode);
  this->update();
}

template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::Synchronize(void) 
{
  //float iwDiff  = this->iwMax() - this->iwMin();
  //float b       = (float)((int)log10(iwDiff)-2);
  //double iwMin  = ((int)(this->iwMin()*pow((float)10, (float)-b)))/pow((float)10,(float)-b);
  //double iwMax  = ((int)(this->iwMax()*pow((float)10, (float)-b)))/pow((float)10,(float)-b);
  //double iwStep = (iwMax-iwMin)/100.0;
  //sliceNumberSlider->range( 0.0f, glSliceView->numSlices() );
  //intensityWindowingMinSlider->range(iwMin-iwStep,iwMax+iwStep);
  //intensityWindowingMaxSlider->range(iwMin-iwStep,iwMax+iwStep);
  //sliceNumberSlider->value((float)glSliceView->sliceNum());
  //intensityWindowingMinSlider->step(iwStep);
  //intensityWindowingMaxSlider->step(iwStep);
  //intensityWindowingMinSlider->value(iwMin);
  //intensityWindowingMaxSlider->value(iwMax);
}



template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::SelectSlice(unsigned int num)
{
  //sliceNumberSlider->value(num);
  //glSliceView->sliceNum((int)sliceNumberSlider->value());
  this->sliceNum(num);
  this->update();
}



template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::SetIntensityWindowingMin(float val)
{
  //intensityWindowingMinSlider->value(val);
  //glSliceView->iwMin(intensityWindowingMinSlider->value());  
  this->iwMin(val);  	
  this->update();
}


template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::SetIntensityWindowingMax(float val)
{
  //intensityWindowingMaxSlider->value(val);
  //glSliceView->iwMax(intensityWindowingMaxSlider->value());  
  this->iwMax(val);  
  this->update();
}


template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::CenterWindow(void)
{
  this->winCenter();
  this->update();
}


template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::ZoomIn(void)
{
  float z = this->winZoom();
  this->winZoom(z*2.0f);
  this->update();
}

template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::ZoomOut(void)
{
  float z = this->winZoom();
  this->winZoom(z/2.0f);
  this->update();
}


template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::ShiftUp(void)
{
  this->winShift(1,0);
  this->update();
}


template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::ShiftDown(void)
{
  this->winShift(-1,0);
  this->update();
}


template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::ShiftLeft(void)
{
  this->winShift(0,-1);
  this->update();
}


template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::ShiftRight(void)
{
  this->winShift(0,1);
  this->update();
}

template <class ImagePixelType, class OverlayPixelType>
void
MyImageViewer<ImagePixelType,OverlayPixelType>
::SetOrientation(int val)
{
  //glSliceView->orientation( orientationChoice->value() );
  this->orientation( val );
  this->update();
  Synchronize();
}


template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::ClickSelectCallBack( void (*newClickSelectArgCallBack)(float, float,
                                                         float, float,int,
                                                         void *),
                                                         void * newClickSelectArg)
{
  this->clickSelectCallBack( newClickSelectArgCallBack, 
                                    newClickSelectArg           ); 
}


template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::ViewDetails(bool detail)
{
  this->viewDetails(detail);
}
  
  
template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::ViewValue(bool value)
{
  this->viewValue(value);
}
  
template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::ViewCrosshairs(bool crosshairs)
{
  this->viewCrosshairs(crosshairs);
}


template <class ImagePixelType, class OverlayPixelType>
unsigned int 
MyImageViewer<ImagePixelType,OverlayPixelType>
::GetNumSlices()
{
  return this->GetSliceNum();
}

template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::ChangeActiveClass(int value)
{
  this->change_activeclass(value);
}

template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::ChangeColorMode(int value)
{
	this->ChangeColorTable(value);
}
	
template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::SetColorMax(double value)
{
	this->setColorMax(value);
	this->update();
}
	
template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::SetColorMin(double value)
{
	this->setColorMin(value);
	this->update();
}
	
	
template <class ImagePixelType, class OverlayPixelType>
void 
MyImageViewer<ImagePixelType,OverlayPixelType>
::ClearClickedPoints(void)
{
  this->clearClickedPointsStored();
}

template <class ImagePixelType, class OverlayPixelType>
int 
MyImageViewer<ImagePixelType,OverlayPixelType>
::handle(int event)
{
  switch(event) {        
    case FL_SHORTCUT:
      if (Fl::event_key() == 'b') {
        //std::cout << "pressed b " << std::endl; 
        this->ClearClickedPoints();
      }
      if (Fl::event_key() == 'n') {
        std::cout << "pressed n" << std::endl; 
        this->deleteLastClickedPointsStored();
      }
      break;      
    default:
      break;	
   }

  return GLSliceView<ImagePixelType, OverlayPixelType>::handle(event);
}


} // end namespace itk

#endif