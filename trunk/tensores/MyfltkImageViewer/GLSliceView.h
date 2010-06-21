/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: GLSliceView.h,v $
 Language:  C++
 Date:      $Date: 2005/12/08 18:20:44 $
 Version:   $Revision: 1.25 $
 
  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
  
   This software is distributed WITHOUT ANY WARRANTY; without even 
   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
   PURPOSE.  See the above copyright notices for more information.
   
=========================================================================*/
#ifndef GLSLICEVIEW_H
#define GLSLICEVIEW_H

#include <FL/gl.h>
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Gl_Window.H>

#include "itkColorTable.h"

#include "SliceView.h"

#include <math.h>

namespace itk {
  
/**
* GLSliceView : Derived from abstract class SliceView and Fl_Gl_Window
* See SliceView.h for details...
  **/
  template <class ImagePixelType, class OverlayPixelType>
    class GLSliceView : public SliceView<ImagePixelType>, 
    public Fl_Gl_Window 
    {
public:
  
  typedef Image<ImagePixelType,3>     ImageType;
  typedef Image<OverlayPixelType,3>   OverlayType;
  typedef typename ImageType::Pointer      ImagePointer;
  typedef typename OverlayType::Pointer    OverlayPointer;
  typedef typename ImageType::RegionType   RegionType;
  typedef typename ImageType::SizeType     SizeType;
  typedef typename ImageType::IndexType    IndexType;
  
  typedef itk::ColorTable<float>                ColorTableType;
  typedef typename ColorTableType::Pointer      ColorTablePointer;
  
  
protected:
  bool        cValidOverlayData;
  float       cOverlayOpacity;
  
  OverlayPointer cOverlayData;
  void     (* cViewOverlayCallBack)(void);
  
  unsigned char * cWinOverlayData;
  

  ColorTablePointer      cColorTable;
  ColorTablePointer      cColorTableOv;		
  ColorTablePointer      cnewColorTable;
  unsigned int           cOverlayColorIndex;
  
public:
/*! FLTK required constructor - must use imData() to complete 
  definition */
  GLSliceView(int x, int y, int w, int h, const char *l);
  
  /*! Specify the 3D image to view slice by slice */
  virtual void SetInputImage(ImageType * newImData);
  virtual const ImagePointer & GetInputImage(void) const;
  
  /*! Specify the 3D image to view as an overlay */
  void SetInputOverlay(OverlayType * newOverlayData);
  
  /*! Return a pointer to the overlay data */
  const OverlayPointer & GetInputOverlay(void) const;
  
  /*! Turn on/off the viewing of the overlay */
  void  ViewOverlayData(bool newViewOverlayData);
  
  /*! Status of the overlay - viewed /not viewed */
  bool  ViewOverlayData(void);
  
  /*! Specify the opacity of the overlay */
  void  OverlayOpacity(float newOverlayOpacity);
  
  /*! Get the opacity of the overlay */
  float OverlayOpacity(void);
  
  /*! Called when overlay is toggled or opacity is changed */
  void  ViewOverlayCallBack(void (* newOverlayCallBack)(void));
  
  void SetBuble(float x, float y, float z, int Radius, float value);
  void ClearOverlay();
		
  ColorTablePointer GetColorTable(void);
  void SetColorTable(ColorTablePointer newColorTable);
  void SetROISColorTable(ColorTablePointer colorTable);
  void SetFAColorTable(ColorTablePointer colorTable);
  void ChangeColorTable( int value );
  void setColorMax( double value ) {
	  this->cColorDataMax = value; 
  };
  void setColorMin( double value ) {
	  this->cColorDataMin = value; 
  };
		
  /*! Turn on/off the display of clicked points */
  void ViewClickedPoints( bool newViewClickedPoints );
 
  unsigned int GetSliceNum( );

  /*! Status of clicked points display - on/off */
  bool ViewClickedPoints();

  void DrawLateralLetters(float* v);
  void DrawVerticalLetters(float* v);

  virtual void clickSelect(float x, float y, float z);
  virtual void clickMove(float x, float y, float z);

  virtual void size(int w, int h);
  virtual void resize(int x, int y, int w, int h);
  
  virtual void update();
  virtual void draw();
  
  virtual int  handle(int event);
  
  /*! Display Overlay in Color 'c'. You must ensure that the color-table specified
   * contains color 'c'. For example with the default useDiscrete() color table,
   * SetOverlayColorIndex( 0 ) will display the overlay in red. 
   * SetOverlayColorIndex( 1 ) purple etc.... */
  void SetOverlayColorIndex( unsigned int c)
    {
    cOverlayColorIndex = c;
    }
  };
  
  
  
template <class ImagePixelType, class OverlayPixelType>
GLSliceView<ImagePixelType, OverlayPixelType>::
GLSliceView(int x, int y, int w, int h, const char *l):
SliceView<ImagePixelType>(x, y, w, h, l), Fl_Gl_Window(x, y, w, h, l)
  {
  when(FL_WHEN_NOT_CHANGED | FL_WHEN_ENTER_KEY);
  cValidOverlayData     = false;
  this->cViewOverlayData      = false;
  cViewOverlayCallBack  = NULL;
  cOverlayOpacity       = 0.0;
  cWinOverlayData       = NULL;
  cColorTable = ColorTableType::New();
  cColorTableOv = ColorTableType::New();
	  
  cOverlayColorIndex = 33;  // este valor es necesario para los colores de overlaydata 
  cColorTableOv->UseHeatColors(33);
  cColorTable->UseHeatColors(62);
  this->SetFAColorTable(cColorTable);
  this->SetROISColorTable(cColorTableOv);
  
  //cnewColorTable = ColorTableType::New();
  //cnewColorTable->useHeat();
  //cnewColorTable->useDiscrete();
  }
  
template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
ChangeColorTable( int value ) {
	switch (value) {
		case 0:
	       cColorTable->UseHeatColors(62);
		   this->SetFAColorTable(cColorTable);
		   break;
		case 1:
		    cColorTable->UseHeatColors(32);
			this->SetROISColorTable(cColorTable);
			break;
		case 2:
			cColorTable->UseHeatColors(62);
			break;
		case 3:
			cColorTable->useDiscrete();
			break;
		default:
			break;
	}
			
}
		
template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
SetROISColorTable( ColorTablePointer colorTable ) {
	colorTable->SetColor(0,1,0,0,"rojo");
	colorTable->SetColor(1,1,1,0.0,"amarillo");
	colorTable->SetColor(2,0.9,0.0,0.5,"fuxia");
	colorTable->SetColor(3,0.2,0.2,0.8,"azul violeta");
	colorTable->SetColor(4,0.1,0.5,0.1,"verde oscuro");
	colorTable->SetColor(5,0,1,0,"verde");
	colorTable->SetColor(6,1,1,0.6,"beige");
	colorTable->SetColor(7,0.9,0.3,0.9,"rosa palido");
	colorTable->SetColor(8,0.7,0.3,0.2,"marron");
	colorTable->SetColor(9,0.9,0.9,0.9,"gris claro");
	colorTable->SetColor(10,0.5,1,0.5,"verde palido");
	colorTable->SetColor(11,0.1,0.4,0.4,"azul verdoso");
	colorTable->SetColor(12,0.6,0.6,0.1,"amarillo palido");
	colorTable->SetColor(13,0,0,1,"azul");
	colorTable->SetColor(14,0.6,0.5,1,"violeta claro");
	colorTable->SetColor(15,0.6,0.1,0.1,"marron oscuro?");
	colorTable->SetColor(16,0.5,0,0,"rojo claro?");
	colorTable->SetColor(17,0.3,0.7,0.3,"verde palido");
	colorTable->SetColor(18,0,1,1,"celeste");
	colorTable->SetColor(19,1,1,1,"blanco");
	colorTable->SetColor(20,0.8,0.2,0.7,"violeta oscuro");
	colorTable->SetColor(21,0.7,0.5,0.2,"marron claro");
	colorTable->SetColor(22,0.2,0.7,0.5,"verde palido");
	colorTable->SetColor(23,0.5,1,1,"azul palido");
	colorTable->SetColor(24,1,0.5,1,"rosa");
	colorTable->SetColor(25,0.5,0,1,"violeta");
	colorTable->SetColor(26,1,0.5,0,"naranja");
	colorTable->SetColor(27,1,0,0.5,"rosa fuerte");
	colorTable->SetColor(28,0,0.5,1,"azul claro");
	colorTable->SetColor(29,0.2,1,0.6,"undef");
	colorTable->SetColor(30,1,0.7,0.1,"undef");
	colorTable->SetColor(31,1,0.8,0.7,"undef");
	colorTable->SetColor(32,0.3,0.6,0.6,"undef");
}

template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
SetFAColorTable( ColorTablePointer colorTable ) {
	
	colorTable->SetColor(0,0  ,       0  ,  0.5625,"undef");
	colorTable->SetColor(1,0  ,       0  ,  0.6250,"undef");
	colorTable->SetColor(2,0  ,       0  ,  0.6875,"undef");
	colorTable->SetColor(3,0  ,       0  ,  0.7500,"undef");
	colorTable->SetColor(4,0  ,       0  ,  0.8125,"undef");
	colorTable->SetColor(5,0  ,       0  ,  0.8750,"undef");
	colorTable->SetColor(6,0  ,       0  ,  0.9375,"undef");
	colorTable->SetColor(7,0  ,       0  ,  1.0000,"undef");
	colorTable->SetColor(8,0  ,  0.0625  ,  1.0000,"undef");
	colorTable->SetColor(9,0  ,  0.1250  ,  1.0000,"undef");
	colorTable->SetColor(10,0  ,  0.1875  ,  1.0000,"undef");
	colorTable->SetColor(11,0  ,  0.2500  ,  1.0000,"undef");
	colorTable->SetColor(12,0  ,  0.3125  ,  1.0000,"undef");
	colorTable->SetColor(13,0  ,  0.3750  ,  1.0000,"undef");
	colorTable->SetColor(14,0  ,  0.4375  ,  1.0000,"undef");
	colorTable->SetColor(15,0  ,  0.5000  ,  1.0000,"undef");
	colorTable->SetColor(16,0  ,  0.5625  ,  1.0000,"undef");
	colorTable->SetColor(17,0  ,  0.6250  ,  1.0000,"undef");
	colorTable->SetColor(17,0  ,  0.6875  ,  1.0000,"undef");
	colorTable->SetColor(18,0  ,  0.7500  ,  1.0000,"undef");
	colorTable->SetColor(19,0  ,  0.8125  ,  1.0000,"undef");
	colorTable->SetColor(20,0  ,  0.8750  ,  1.0000,"undef");
	colorTable->SetColor(21,0  ,  0.9375  ,  1.0000,"undef");
	colorTable->SetColor(22,0  ,  1.0000  ,  1.0000,"undef");
    colorTable->SetColor(23,0.0625  ,  1.0000  ,  0.9375,"undef");
    colorTable->SetColor(24,0.1250  ,  1.0000  ,  0.8750,"undef");
    colorTable->SetColor(25,0.1875  ,  1.0000  ,  0.8125,"undef");
    colorTable->SetColor(26,0.2500  ,  1.0000  ,  0.7500,"undef");
    colorTable->SetColor(27,0.3125  ,  1.0000  ,  0.6875,"undef");
    colorTable->SetColor(28,0.3750  ,  1.0000  ,  0.6250,"undef");
    colorTable->SetColor(29,0.4375  ,  1.0000  ,  0.5625,"undef");
    colorTable->SetColor(30,0.5000  ,  1.0000  ,  0.5000,"undef");
    colorTable->SetColor(31,0.5625  ,  1.0000  ,  0.4375,"undef");
    colorTable->SetColor(32,0.6250  ,  1.0000  ,  0.3750,"undef");
    colorTable->SetColor(33,0.6875  ,  1.0000  ,  0.3125,"undef");
    colorTable->SetColor(34,0.7500  ,  1.0000  ,  0.2500,"undef");
    colorTable->SetColor(35,0.8125  ,  1.0000  ,  0.1875,"undef");
    colorTable->SetColor(36,0.8750  ,  1.0000  ,  0.1250,"undef");
    colorTable->SetColor(37,0.9375  ,  1.0000  ,  0.0625,"undef");
    colorTable->SetColor(38,1.0000  ,  1.0000  ,       0,"undef");
    colorTable->SetColor(39,1.0000  ,  0.9375  ,       0,"undef");
    colorTable->SetColor(40,1.0000  ,  0.8750  ,       0,"undef");
    colorTable->SetColor(41,1.0000  ,  0.8125  ,       0,"undef");
    colorTable->SetColor(42,1.0000  ,  0.7500  ,       0,"undef");
    colorTable->SetColor(43,1.0000  ,  0.6875  ,       0,"undef");
    colorTable->SetColor(44,1.0000  ,  0.6250  ,       0,"undef");
    colorTable->SetColor(45,1.0000  ,  0.5625  ,       0,"undef");
    colorTable->SetColor(46,1.0000  ,  0.5000  ,       0,"undef");
    colorTable->SetColor(47,1.0000  ,  0.4375  ,       0,"undef");
    colorTable->SetColor(48,1.0000  ,  0.3750  ,       0,"undef");
    colorTable->SetColor(49,1.0000  ,  0.3125  ,       0,"undef");
    colorTable->SetColor(50,1.0000  ,  0.2500  ,       0,"undef");
    colorTable->SetColor(51,1.0000  ,  0.1875  ,       0,"undef");
    colorTable->SetColor(52,1.0000  ,  0.1250  ,       0,"undef");
    colorTable->SetColor(53,1.0000  ,  0.0625  ,       0,"undef");
    colorTable->SetColor(54,1.0000  ,       0  ,       0,"undef");
    colorTable->SetColor(55,0.9375  ,       0  ,       0,"undef");
    colorTable->SetColor(56,0.8750  ,       0  ,       0,"undef");
    colorTable->SetColor(57,0.8125  ,       0  ,       0,"undef");
    colorTable->SetColor(58,0.7500  ,       0  ,       0,"undef");
    colorTable->SetColor(59,0.6875  ,       0  ,       0,"undef");
    colorTable->SetColor(60,0.6250  ,       0  ,       0,"undef");
    colorTable->SetColor(61,0.5625  ,       0  ,       0,"undef");
    colorTable->SetColor(62,0.5000  ,       0  ,       0,"undef");
}		
	
//
// Set the input image to be displayed
// Warning: the current overlay is destroyed if the size of the image
// is different from the size of the overlay.
//
template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
SetInputImage(ImageType * newImData)
  {
  RegionType region = newImData->GetLargestPossibleRegion();
  if( region.GetNumberOfPixels() == 0 ) 
    {
    return;
    }

  SizeType   imSize   = region.GetSize();

  // If the overlay has been set and the size is different from the new image,
  // it is removed.
  if( cValidOverlayData)
    {  
    SizeType overlay_size = cOverlayData->GetLargestPossibleRegion().GetSize();
    
    if((overlay_size[0] != imSize[0])
      ||  (overlay_size[1] != imSize[1])
      ||  (overlay_size[2] != imSize[2]))
      {
       if(cWinOverlayData != NULL)
         {
         delete [] cWinOverlayData;
         }
       cWinOverlayData       = NULL;
       cValidOverlayData     = false;
      }       
    }

  this->cImData = newImData;
  this->cDimSize[0]=imSize[0];
  this->cDimSize[1]=imSize[1];
  this->cDimSize[2]=imSize[2];
  this->cSpacing[0]=this->cImData->GetSpacing()[0];
  this->cSpacing[1]=this->cImData->GetSpacing()[1];
  this->cSpacing[2]=this->cImData->GetSpacing()[2];
  this->cOrigin[0]=this->cImData->GetOrigin()[0];
  this->cOrigin[1]=this->cImData->GetOrigin()[1];
  this->cOrigin[2]=this->cImData->GetOrigin()[2];
    
  //calculating cDataMax and cDataMin    
  IndexType ind;
  ind[0] = 0; 
  ind[1] = 0; 
  ind[2] = 0;
  
  this->cDataMax = this->cImData->GetPixel(ind);
  this->cDataMin = this->cDataMax;
  ImagePixelType tf;
  
  
  for( unsigned int i=0; i<this->cDimSize[0]; i++ )
    {
    ind[0] = i;
    for(unsigned int j=0; j<this->cDimSize[1]; j++ )
      {
      ind[1] = j;
      for( unsigned int k=0; k<this->cDimSize[2]; k++ )
        {
        ind[2] = k;
        tf = this->cImData->GetPixel(ind);
        if(tf > this->cDataMax) 
          {
          this->cDataMax = tf;
          }
        else 
          {
          if(tf < this->cDataMin)
            {
            this->cDataMin = tf;
            }
          }
        }
      }
    }
  
  this->cIWMin      = this->cDataMin;
  this->cIWMax      = this->cDataMax;
  this->cIWModeMin  = IW_MIN;
  this->cIWModeMax  = IW_MAX;
  this->cColorDataMax = this->cDataMax;
  this->cColorDataMin = this->cDataMin;
	  
  this->cImageMode = IMG_VAL;
  
  this->cWinZoom = 1;
  
  this->cWinOrientation = 2;
  this->cWinOrder[0] = 0;
  this->cWinOrder[1] = 1;
  this->cWinOrder[2] = 2;
  
  this->cWinCenter[0] = this->cDimSize[0]/2;
  this->cWinCenter[1] = this->cDimSize[1]/2;
  this->cWinCenter[2] = 0;
  
  this->cWinMinX  = 0;
  this->cWinSizeX = this->cDimSize[0];
  if(this->cWinSizeX<this->cDimSize[1])
    {
    this->cWinSizeX = this->cDimSize[1];
    }
  if(this->cWinSizeX<this->cDimSize[2])
    {
    this->cWinSizeX = this->cDimSize[2];
    }
  this->cWinMaxX  = this->cWinSizeX - 1;
  
  this->cWinMinY  = 0;
  this->cWinSizeY = this->cWinSizeX;
  this->cWinMaxY  = this->cWinSizeY - 1;
  
  this->cWinDataSizeX = this->cWinMaxX;//this->cDimSize[0];
  this->cWinDataSizeY = this->cWinMaxY;//this->cDimSize[1];
  
  if(this->cWinImData != NULL)
    {
    delete [] this->cWinImData;
    }
  this->cWinImData = new unsigned char[ this->cWinDataSizeX * this->cWinDataSizeY ];
  
  if(this->cWinImColorData != NULL)
    {
    delete [] this->cWinImColorData;
    }
  this->cWinImColorData = new unsigned char[ this->cWinDataSizeX * this->cWinDataSizeY * 3];

  if(this->cWinZBuffer != NULL) 
    {
    delete [] this->cWinZBuffer;
    }
  this->cWinZBuffer = new unsigned short[ this->cWinDataSizeX * this->cWinDataSizeY ];
  
  this->cViewImData  = true;
  this->cValidImData = true;
  
  }





//
//
template <class ImagePixelType, class OverlayPixelType>
const typename Image<ImagePixelType,3>::Pointer &
GLSliceView<ImagePixelType, OverlayPixelType>
::GetInputImage(void) const
  {
  return this->cImData;
  }



//
//
template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>
::SetInputOverlay( OverlayType * newOverlayData )
  {
  RegionType newoverlay_region = 
    newOverlayData->GetLargestPossibleRegion();

  SizeType   newoverlay_size   = newoverlay_region.GetSize();

  if( this->cValidImData 
      &&  (newoverlay_size[0] == this->cDimSize[0])
      &&  (newoverlay_size[1] == this->cDimSize[1])
      &&  (newoverlay_size[2] == this->cDimSize[2])
    )
    {
    this->cOverlayData = newOverlayData;
    
    this->cViewOverlayData  = true;
    this->cValidOverlayData = true;
    this->cOverlayOpacity   = (float)1.0;
    
    if(this->cWinOverlayData != NULL) 
      {
      delete [] this->cWinOverlayData;
      }
    

    const unsigned long bufferSize = this->cWinDataSizeX * this->cWinDataSizeY * 4;
    this->cWinOverlayData = new unsigned char[ bufferSize ];
    }
  else // return a warning
    {
      if(!this->cValidImData)
        {
        std::cout << "GLSliceView: Please set the input image before overlay"  
                  << std::endl;
        std::cout << "GLSliceView: Overlay not set." << std::endl;
        }
      else if((newoverlay_size[0] != this->cDimSize[0])
      ||  (newoverlay_size[1] != this->cDimSize[1])
      ||  (newoverlay_size[2] != this->cDimSize[2])
      )
        {
        std::cout << "GLSliceView: Error: overlay and input images should be the same size" 
                  << std::endl;
        std::cout << "GLSliceView: Overlay not set." << std::endl;
        }
      
    }
  }

template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
SetBuble(float x, float y, float z, int Radius, float value) {
	if (this->cValidOverlayData) {
	  IndexType pixelIndex;
	  //std::cout << " Pintando en x " << x << " , y, " << y << " , z, " << z << " con Radius " << Radius << " con valor " << value << std::endl;
	 if (Radius == 0) {
		 pixelIndex[0] = (int)x;
		 pixelIndex[1] = (int)y;
		 pixelIndex[2] = (int)z;
		 cOverlayData->SetPixel(pixelIndex, value);
	 } else {
			
	  for(int i=-Radius; i<= Radius; i++) {
		for(int j=-Radius; j<= Radius; j++) {
		  //for(int k=-Radius; k<= Radius; k++) {
			if (i*i + j*j <= Radius*Radius) {
			  pixelIndex[0] = (int)x+i;
			  pixelIndex[1] = (int)y+j;
			  //pixelIndex[2] = (int)z+k;
			  pixelIndex[2] = (int)z;
			  cOverlayData->SetPixel(pixelIndex, value);
			}
		  //}
	    }
       }
     }
	}
	this->update();
}
	
template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
ClearOverlay( ) {
  cOverlayData->FillBuffer(0);
}

template <class ImagePixelType, class OverlayPixelType>
const typename GLSliceView<ImagePixelType, 
OverlayPixelType>::OverlayPointer &
GLSliceView<ImagePixelType, OverlayPixelType>::GetInputOverlay( void ) 
const
  {
  return this->cOverlayData;
  }




//
//
template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
ViewOverlayData( bool newViewOverlayData)
  {
  
  this->cViewOverlayData = newViewOverlayData;
  
  if( this->cViewOverlayCallBack != NULL )
    {
    cViewOverlayCallBack();
    }
  
  this->redraw();
  
  }

template <class ImagePixelType, class OverlayPixelType>
bool 
GLSliceView<ImagePixelType, OverlayPixelType>::
ViewOverlayData(void)
  {
  
  return this->cViewOverlayData;
  
  }


template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
ViewOverlayCallBack( void (* newViewOverlayCallBack)(void) )
  {
  this->cViewOverlayCallBack = newViewOverlayCallBack;
  }


template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
ViewClickedPoints( bool newViewClickedPoints )
{
    this->cViewClickedPoints = newViewClickedPoints;

    this->redraw();
}

template <class ImagePixelType, class OverlayPixelType>
bool
GLSliceView<ImagePixelType, OverlayPixelType>::
ViewClickedPoints()
{
    return this->cViewClickedPoints;
}

template <class ImagePixelType, class OverlayPixelType>
unsigned int
GLSliceView<ImagePixelType, OverlayPixelType>::
GetSliceNum()
{
    return this->sliceNum();
}

//
//
template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
OverlayOpacity( float newOverlayOpacity )
  {
  this->cOverlayOpacity = newOverlayOpacity;
  
  if(this->cViewOverlayCallBack != NULL) 
    {
    this->cViewOverlayCallBack();
    }
  }




template <class ImagePixelType, class OverlayPixelType>
float 
GLSliceView<ImagePixelType, OverlayPixelType>::
OverlayOpacity(void)
  {
  return this->cOverlayOpacity;
  }




//
//
//
template <class ImagePixelType, class OverlayPixelType>
typename GLSliceView<ImagePixelType, OverlayPixelType>::ColorTablePointer 
GLSliceView<ImagePixelType, OverlayPixelType>::
GetColorTable(void)
  {
  return this->cColorTable;
  }



//
//
//
template <class ImagePixelType, class OverlayPixelType>
void
GLSliceView<ImagePixelType, OverlayPixelType>::
SetColorTable(typename 
              GLSliceView<ImagePixelType, OverlayPixelType>::ColorTablePointer 
              newColorTable)
  {
  cColorTable = newColorTable;
  }



//
//
template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
update()
  {
  
  if( !this->cValidImData ) 
    {
    return;
    }

  double zoomBase = this->cW / (this->cDimSize[this->cWinOrder[0]]
                                * (fabs(this->cSpacing[this->cWinOrder[0]])
                                   /fabs(this->cSpacing[0])));
  if(zoomBase > this->cH / (this->cDimSize[this->cWinOrder[1]]
                            * (fabs(this->cSpacing[this->cWinOrder[1]])
                              / fabs(this->cSpacing[0]))))
    {
    zoomBase = this->cH / (this->cDimSize[this->cWinOrder[1]]
                           * (fabs(this->cSpacing[this->cWinOrder[1]])
                              / fabs(this->cSpacing[0])));
    }
  double scale0 = this->cWinZoom * zoomBase 
                                 * fabs(this->cSpacing[this->cWinOrder[0]])
                                 / fabs(this->cSpacing[0]);
  double scale1 = this->cWinZoom * zoomBase 
                                 * fabs(this->cSpacing[this->cWinOrder[1]]) 
                                 / fabs(this->cSpacing[0]);

  if(this->cWinZoom>1)
    {
    this->cWinSizeX = (int)( this->cW / scale0 );
    this->cWinMinX = (int)( (int)this->cWinCenter[ this->cWinOrder[0] ] 
                            - this->cWinSizeX/2 );
    this->cWinMaxX = (int)( (int)this->cWinCenter[ this->cWinOrder[0] ] 
                            + this->cWinSizeX/2 );
    }
  else
    {
    this->cWinSizeX = (int)(this->cDimSize[ this->cWinOrder[0] ]);
    this->cWinMinX = 0;
    this->cWinMaxX = (int)( (int)(this->cDimSize[ this->cWinOrder[0] ]) - 1 );
    this->cWinCenter[this->cWinOrder[0]] = 
                     (int)( this->cDimSize[ this->cWinOrder[0] ] / 2);
    }
  if( this->cWinMinX <= - (int) this->cDimSize[ this->cWinOrder[0] ] ) 
    {
    this->cWinMinX = -(int)this->cDimSize[ this->cWinOrder[0] ] + 1;
    }
  else if(this->cWinMinX >= (int)this->cDimSize[ this->cWinOrder[0] ]) 
    {
    this->cWinMinX = (int)this->cDimSize[ this->cWinOrder[0] ] - 1;
    }
  if( this->cWinMaxX >= (int)( this->cDimSize[ this->cWinOrder[0] ] ) )
    {
    this->cWinMaxX = (int)this->cDimSize[ this->cWinOrder[0] ] - 1;
    }
  
  if(this->cWinZoom>1)
    {
    this->cWinSizeY = (int)( this->cH / scale1 );
    this->cWinMinY = (int)( (int)(this->cWinCenter[ this->cWinOrder[1] ]) 
                             - this->cWinSizeY/2 );
    this->cWinMaxY = (int)( (int)(this->cWinCenter[ this->cWinOrder[1] ]) 
                             + this->cWinSizeY/2 );
    }
  else
    {
    this->cWinSizeY = (int)(this->cDimSize[ this->cWinOrder[1] ]);
    this->cWinMinY = 0;
    this->cWinMaxY = (int)( (int)(this->cDimSize[ this->cWinOrder[1] ]) - 1 );
    this->cWinCenter[this->cWinOrder[1]] = 
                     (int)( this->cDimSize[ this->cWinOrder[1] ] / 2);
    }
  if( this->cWinMinY <= - (int)( this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    this->cWinMinY = -(int)this->cDimSize[ this->cWinOrder[1] ] + 1;
    }
  else if( this->cWinMinY >= (int)(this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    this->cWinMinY = this->cDimSize[ this->cWinOrder[1] ] - 1;
    } 
  if( this->cWinMaxY >= (int)( this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    this->cWinMaxY = this->cDimSize[ this->cWinOrder[1] ] - 1;
    }
  
  memset( this->cWinImData, 0, this->cWinDataSizeX*this->cWinDataSizeY );
  memset( this->cWinImColorData, 0, this->cWinDataSizeX*this->cWinDataSizeY*3 );
  if( this->cValidOverlayData ) 
    {
    memset(this->cWinOverlayData, 0, this->cWinDataSizeX*this->cWinDataSizeY*4);
    }
  
  IndexType ind;
  
  int l, m;
  
  float tf;
  
  ind[ this->cWinOrder[ 2 ] ] = this->cWinCenter[ this->cWinOrder[ 2 ] ];
  int startK = this->cWinMinY;
  if(startK<0)
    startK = 0;
  int startJ = this->cWinMinX;
  if(startJ<0)
    startJ = 0;
  for(int k=startK; k <= this->cWinMaxY; k++)
    {
    ind[this->cWinOrder[1]] = k;
    
    if(k-this->cWinMinY >= (int)this->cWinDataSizeY)
      continue;

    for(int j=startJ; j <= this->cWinMaxX; j++) 
      {
      ind[this->cWinOrder[0]] = j;
      
      if(j-this->cWinMinX >= (int)this->cWinDataSizeX)
         continue;

      switch( this->cImageMode ) 
        {
        default:
        case IMG_VAL:
          tf = (float)((this->cImData->GetPixel(ind)-this->cIWMin) 
                       / (this->cIWMax-this->cIWMin)*255);
          break;
        case IMG_INV:
          tf = (float)((this->cIWMax-this->cImData->GetPixel(ind)) 
                       / (this->cIWMax-this->cIWMin)*255);
          break;
        case IMG_LOG:
          tf = (float)(log(this->cImData->GetPixel(ind)-this->cIWMin+0.00000001)
                       /log(this->cIWMax-this->cIWMin+0.00000001)*255);
          break;
        case IMG_DX:
          if(ind[0]>0) 
            {
            tf = (float)((this->cImData->GetPixel(ind)-this->cIWMin)
                         / (this->cIWMax-this->cIWMin)*255);
            ind[0]--;
            tf -= (float)((this->cImData->GetPixel(ind)-this->cIWMin) 
                         / (this->cIWMax-this->cIWMin)*255);
            ind[0]++;
            tf += 128;
            } 
          else
            {
            tf = 128;
            }
          break;
        case IMG_DY:
          if(ind[1]>0) 
            {
            tf = (float)((this->cImData->GetPixel(ind)-this->cIWMin) 
                         / (this->cIWMax-this->cIWMin)*255);
            ind[1]--;
            tf -= (float)((this->cImData->GetPixel(ind)-this->cIWMin) 
                          / (this->cIWMax-this->cIWMin)*255);
            ind[1]++;
            tf += 128;
            }
          else
            {
            tf = 128;
            }
          break;
        case IMG_DZ:
          if(ind[2]>0) 
            {
            tf = (float)((this->cImData->GetPixel(ind)-this->cIWMin)
                         / (this->cIWMax-this->cIWMin)*255);
            ind[2]--;
            tf -= (float)((this->cImData->GetPixel(ind)-this->cIWMin)
                          / (this->cIWMax-this->cIWMin)*255);
            ind[2]++;
            tf += 128;
            }
          else
            {
            tf = 128;
            }
          break;
        case IMG_BLEND:
          {
          const int tempval = (int)this->cWinCenter[this->cWinOrder[2]]-1;
          int tmpI = ind[this->cWinOrder[2]];
          ind[this->cWinOrder[2]] = (tempval < 0 ) ? 0 : tempval;
          tf = (float)(this->cImData->GetPixel(ind));
          
          ind[this->cWinOrder[2]] = this->cWinCenter[this->cWinOrder[2]];
          tf += (float)(this->cImData->GetPixel(ind))*2;
          
          const int tempval1 = (int)this->cDimSize[this->cWinOrder[2]]-1;
          const int tempval2 = (int)this->cWinCenter[this->cWinOrder[2]]+1;
          ind[this->cWinOrder[2]] = (tempval1 < tempval2 ) ? tempval1 : tempval2;
          tf += (float)(this->cImData->GetPixel(ind));
          
          tf = (float)((tf/4-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          ind[this->cWinOrder[2]] = tmpI;
          break;
          }
        case IMG_MIP:
          {
          tf = this->cIWMin;
          m = (j-this->cWinMinX) + (k-this->cWinMinY)*this->cWinDataSizeX;
          this->cWinZBuffer[m] = 0;
          int tmpI = ind[this->cWinOrder[2]];
          for(l=0; l<(int)this->cDimSize[this->cWinOrder[2]]; l++) 
            {
            ind[this->cWinOrder[2]] = l;        
            if(this->cImData->GetPixel(ind) > tf) 
              {
              tf = (float)(this->cImData->GetPixel(ind));
              this->cWinZBuffer[m] = (unsigned short)l;
              }
            }
          tf = (float)((tf-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          ind[this->cWinOrder[2]] = tmpI;
          break;       
	  } 
        case IMG_COLOR:
          {
	  tf = (float)(this->cImData->GetPixel(ind));	  
	  break;
	  }
	}

      
      if( tf > 255 )
        {
        switch(this->cIWModeMax) 
          {
          case IW_MIN:
            tf = 0;
            break;
          default:
          case IW_MAX:
            tf = 255;
            break;
          case IW_FLIP:
            tf = 512-tf;
            if(tf<0) 
              {
              tf = 0;
              }
            break;
          }
        }
      else 
        {
        if( tf < 0 )
          {
          switch(this->cIWModeMin) 
            {
            default:
            case IW_MIN:
              tf = 0;
              break;
            case IW_MAX:
              tf = 255;
              break;
            case IW_FLIP:
              tf = -tf;
              if(tf>255)
                {
                tf = 255;
                }
              break;
            }
          }
        }
      
      l = (j-this->cWinMinX) + (k-this->cWinMinY)*this->cWinDataSizeX;
      this->cWinImData[l] = (unsigned char)tf;
      unsigned int overlayColorIndex = 0;
		  float f;
      if (this->cViewColorMode ) {
	    int n = l*3;
	    if (sizeof(ImagePixelType) == 1) {
            m = (int)*((unsigned char *)&(this->cImData->GetPixel(ind)));
		} 
		if (sizeof(ImagePixelType) == 2) { 
            m = (int)*((unsigned short *)&(this->cImData->GetPixel(ind)));
		} 
		if (sizeof(ImagePixelType) == 4) { 
            //m = (int)*((float *)&(this->cImData->GetPixel(ind)));
			f= this->cImData->GetPixel(ind)*cColorTable->GetNumberOfColors()/(this->cColorDataMax-this->cColorDataMin);
			//f= this->cImData->GetPixel(ind)*cColorTable->GetNumberOfColors()/(this->cDataMax-this->cDataMin);
			//f= this->cImData->GetPixel(ind)*cColorTable->GetNumberOfColors()/this->cDataMax;
			m = (int)*((float *)&(f));
		}
		if( m >= (int)cColorTable->GetNumberOfColors() )  { 
			m = cColorTable->GetNumberOfColors();
		  //m = m % cColorTable->GetNumberOfColors();
	    }
		  
		if( m > 0 ) {
            overlayColorIndex = m-1;
            this->cWinImColorData[n+0] = 
              (unsigned char)(cColorTable->GetColorComponent(overlayColorIndex,
															 'r') * 255);
		    this->cWinImColorData[n+1] = 
              (unsigned char)(cColorTable->GetColorComponent(overlayColorIndex,
                                                             'g') * 255);
		    this->cWinImColorData[n+2] = 
              (unsigned char)(cColorTable->GetColorComponent(overlayColorIndex,
                                                             'b') * 255);
			
		}
      }

      if( this->cValidOverlayData ) {

        l = l * 4;
        if(this->cImageMode == IMG_MIP) {
          ind[this->cWinOrder[2]] = this->cWinZBuffer[(j-this->cWinMinX) + 
			           (k-this->cWinMinY)*this->cWinDataSizeX];
	    }
        
        if( sizeof( OverlayPixelType ) == 1  ||
            sizeof( OverlayPixelType ) == 2      )
          {
          if (sizeof( OverlayPixelType ) == 1)
            {
            m = (int)*((unsigned char *)&(cOverlayData->GetPixel(ind)));
            }
          else
            {
            m = (int)*((unsigned short *)&(cOverlayData->GetPixel(ind)));
            }
          if( m >= (int)cColorTableOv->GetNumberOfColors() ) 
            { 
            m = cColorTableOv->GetNumberOfColors() - 1;
            }
          if( m > 0 ) {
            overlayColorIndex = m-1;
			if( static_cast<unsigned int>(m) > cOverlayColorIndex )
              {
              overlayColorIndex = cOverlayColorIndex;
              }
            cWinOverlayData[l+0] = 
              (unsigned char)(cColorTableOv->GetColorComponent(overlayColorIndex,
                                                             'r') * 255);
            cWinOverlayData[l+1] = 
              (unsigned char)(cColorTableOv->GetColorComponent(overlayColorIndex,
                                                             'g') * 255);
            cWinOverlayData[l+2] = 
              (unsigned char)(cColorTableOv->GetColorComponent(overlayColorIndex,
                                                             'b') * 255);
            cWinOverlayData[l+3] = 
              (unsigned char)(cOverlayOpacity*255);
            }
          }
        else 
          {
          if(((unsigned char *)&(cOverlayData->GetPixel(ind)))[0]
            + ((unsigned char *)&(cOverlayData->GetPixel(ind)))[1]
            + ((unsigned char *)&(cOverlayData->GetPixel(ind)))[2] > 0)
            {
            if( sizeof( OverlayPixelType ) == 3 )
              {
              cWinOverlayData[l+0] = 
                ((unsigned char *)&(cOverlayData->GetPixel(ind)))[0];
              cWinOverlayData[l+1] = 
                ((unsigned char *)&(cOverlayData->GetPixel(ind)))[1];
              cWinOverlayData[l+2] = 
                ((unsigned char *)&(cOverlayData->GetPixel(ind)))[2];
              cWinOverlayData[l+3] = 
                (unsigned char)(cOverlayOpacity*255);
              }
            else 
              {
              if( sizeof( OverlayPixelType ) == 4 ) 
                {
                cWinOverlayData[l+0] = 
                  ((unsigned char *)&(cOverlayData->GetPixel(ind)))[0];
                cWinOverlayData[l+1] = 
                  ((unsigned char *)&(cOverlayData->GetPixel(ind)))[1];
                cWinOverlayData[l+2] = 
                  ((unsigned char *)&(cOverlayData->GetPixel(ind)))[2];
                cWinOverlayData[l+3] = 
                  (unsigned char)(((unsigned char *)
                  &(cOverlayData->GetPixel(ind)))[3]*cOverlayOpacity);
                }
              }
            }
          }
      } /* endif this->cValidOverlayData */
    }
  }
  
  this->redraw();
  
  
}




template <class ImagePixelType, class OverlayPixelType>
void GLSliceView<ImagePixelType, OverlayPixelType>::
clickSelect(float x, float y, float z)
  {
  SliceView<ImagePixelType>::clickSelect(x, y, z);
  if(this->cViewValue || this->cViewCrosshairs)
    {
    this->redraw();
    }
  }

template <class ImagePixelType, class OverlayPixelType>
void GLSliceView<ImagePixelType, OverlayPixelType>::
clickMove(float x, float y, float z)
  {
    SliceView<ImagePixelType>::clickMove(x ,y ,z);
    if(this->cViewValue || this->cViewCrosshairs)
    {
      this->redraw();
    }
  }

template <class ImagePixelType, class OverlayPixelType>
void GLSliceView<ImagePixelType, OverlayPixelType>::size(int w, int h)
  {
  SliceView<ImagePixelType>::size(w, h);
  Fl_Gl_Window::size(w, h);
  this->update();
  this->redraw();
  }


template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
resize(int x, int y, int w, int h)
  {
  SliceView<ImagePixelType>::resize(x, y, w, h);
  Fl_Gl_Window::resize(x, y, w, h);
  this->update();
  this->redraw();
  }


template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
DrawLateralLetters(float* v)
  {
  
  char sL[2] = "L";
  char sR[2] = "R";
  char sA[2] = "A";
  char sP[2] = "P";
  char sS[2] = "S";
  char sI[2] = "I";

   
  if (fabs(v[0]) > fabs(v[1]) && fabs(v[0]) > fabs(v[2])) {
      if ( (v[0] > 0 && this->cFlipX[this->cWinOrientation] == false) 
	       || (v[0] < 0 && this->cFlipX[this->cWinOrientation] == true) ) { // R L
	    gl_draw( sR, (int)(gl_width(sR)-6), this->cH/2);
        gl_draw( sL, (int)(this->cW-(gl_width(sL)+2)), this->cH/2);
	  } else { // L R
		gl_draw( sL, (int)(gl_width(sL)-6), this->cH/2);
        gl_draw( sR, (int)(this->cW-(gl_width(sR)+2)), this->cH/2);
	  }
   }   
    
	if (fabs(v[1]) > fabs(v[0]) && fabs(v[1]) > fabs(v[2])) {
      if ( (v[1] > 0 && this->cFlipX[this->cWinOrientation] == false) 
	       || (v[1] < 0 && this->cFlipX[this->cWinOrientation] == true) ) { // A P
	    gl_draw( sA, (int)(gl_width(sA)-6), this->cH/2);
        gl_draw( sP, (int)(this->cW-(gl_width(sP)+2)), this->cH/2);
	  } else { // P A
		gl_draw( sP, (int)(gl_width(sP)-6), this->cH/2);
        gl_draw( sA, (int)(this->cW-(gl_width(sA)+2)), this->cH/2);
	  }
   }   
   
   if (fabs(v[2]) > fabs(v[0]) && fabs(v[2]) > fabs(v[1])) {
      if ( (v[2] > 0 && this->cFlipX[this->cWinOrientation] == false) 
	       || (v[2] < 0 && this->cFlipX[this->cWinOrientation] == true) ) { // I S
	    gl_draw( sI, (int)(gl_width(sI)-6), this->cH/2);
        gl_draw( sS, (int)(this->cW-(gl_width(sS)+2)), this->cH/2);
	  } else { // S I
		gl_draw( sS, (int)(gl_width(sS)-6), this->cH/2);
        gl_draw( sI, (int)(this->cW-(gl_width(sI)+2)), this->cH/2);
	  }
   }       
		
}

template <class ImagePixelType, class OverlayPixelType>
void 
GLSliceView<ImagePixelType, OverlayPixelType>::
DrawVerticalLetters(float* v)
  {
  char sL[2] = "L";
  char sR[2] = "R";
  char sA[2] = "A";
  char sP[2] = "P";
  char sS[2] = "S";
  char sI[2] = "I";

   if (fabs(v[0]) > fabs(v[1]) && fabs(v[0]) > fabs(v[2])) {
      if ( (v[0] > 0 && this->cFlipY[this->cWinOrientation] == false) 
	       || (v[0] < 0 && this->cFlipY[this->cWinOrientation] == true)  ) { // L R
	    gl_draw( sL, (int)(this->cW/2), this->cH-12);
	    gl_draw( sR, (int)(this->cW/2), 2 );
	  } else { // R L
	    gl_draw( sR, (int)(this->cW/2), this->cH-12);
	    gl_draw( sL, (int)(this->cW/2), 2 );
	  }
   }   

   if (fabs(v[1]) > fabs(v[0]) && fabs(v[1]) > fabs(v[2])) {
      if ( (v[1] > 0 && this->cFlipY[this->cWinOrientation] == false) 
	       || (v[1] < 0 && this->cFlipY[this->cWinOrientation] == true) ) { // P A
	    gl_draw( sP, (int)(this->cW/2), this->cH-12);
	    gl_draw( sA, (int)(this->cW/2), 2 );
	  } else { // A P
	    gl_draw( sA, (int)(this->cW/2), this->cH-12);
	    gl_draw( sP, (int)(this->cW/2), 2 );
	  }
   }   
  
  if (fabs(v[2]) > fabs(v[0]) && fabs(v[2]) > fabs(v[1])) {
      if ( (v[2] > 0 && this->cFlipY[this->cWinOrientation] == false) 
	       || (v[2] < 0 && this->cFlipY[this->cWinOrientation] == true) ) { // S I
	    gl_draw( sS, (int)(this->cW/2), this->cH-12);
	    gl_draw( sI, (int)(this->cW/2), 2 );
	  } else { // I S
	    gl_draw( sI, (int)(this->cW/2), this->cH-12);
	    gl_draw( sS, (int)(this->cW/2), 2 );
	  }
   }   

}


template <class ImagePixelType, class OverlayPixelType>
void GLSliceView<ImagePixelType, OverlayPixelType>::draw(void)
  {
  if( !valid() )
    {
    glClearColor((float)0.0, (float)0.0, (float)0.0, (float)0.0);          
    glShadeModel(GL_FLAT);
    
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);  //if you don't include this
                                            //image size differences distort
    //glPixelStorei(GL_PACK_ALIGNMENT, 1);
    }
  else
    {
    glClear(GL_COLOR_BUFFER_BIT);    //this clears and paints to black
    
    glMatrixMode(GL_MODELVIEW);    //clear previous 3D draw params
    glLoadIdentity();
    
    glMatrixMode(GL_PROJECTION);
    ortho();
    
    if( !this->cImData ) 
      {
      return;
      }
    
    double zoomBase = this->cW / (this->cDimSize[this->cWinOrder[0]]
                                  * (fabs(this->cSpacing[this->cWinOrder[0]])
                                     / fabs(this->cSpacing[0])));
    if(zoomBase > this->cH / (this->cDimSize[this->cWinOrder[1]]
                              * (fabs(this->cSpacing[this->cWinOrder[1]])
                                 / fabs(this->cSpacing[0]))))
      {
      zoomBase = this->cH / (this->cDimSize[this->cWinOrder[1]]
                             * (fabs(this->cSpacing[this->cWinOrder[1]])
                                / fabs(this->cSpacing[0])));
      }

    double scale0 = this->cWinZoom * zoomBase 
                                   * fabs(this->cSpacing[this->cWinOrder[0]])
                                   / fabs(this->cSpacing[0]);
    double scale1 = this->cWinZoom * zoomBase 
                                   * fabs(this->cSpacing[this->cWinOrder[1]])
                                   / fabs(this->cSpacing[0]);
								   
    int originX = 0;
    int originY = 0;
    if(this->cWinZoom<=1)
      {
      if(this->cW-scale0*this->cDimSize[this->cWinOrder[0]]>0)
        {
        originX = (int)((this->cW-scale0*this->cDimSize[this->cWinOrder[0]])/2.0);
        }
      if(this->cH-scale1*this->cDimSize[this->cWinOrder[1]]>0)
        {
        originY = (int)((this->cH-scale1*this->cDimSize[this->cWinOrder[1]])/2.0);
        }
      }
    glRasterPos2i((this->cFlipX[this->cWinOrientation])?this->cW-originX:originX,
      (this->cFlipY[this->cWinOrientation])?this->cH-originY:originY);  
    glPixelZoom((this->cFlipX[this->cWinOrientation])?-scale0:scale0,
      (this->cFlipY[this->cWinOrientation])?-scale1:scale1);
    
    if( this->cValidImData && this->cViewImData )
      {
		  		  
	  glDrawPixels( this->cWinDataSizeX, this->cWinDataSizeY, 
					   GL_LUMINANCE, GL_UNSIGNED_BYTE, 
					   this->cWinImData );
		 
      }
	  
    if( this->cValidImData && this->cViewImData && this->cViewColorMode)
      {
		
	  glDrawPixels( this->cWinDataSizeX, this->cWinDataSizeY, 
                    GL_RGB, GL_UNSIGNED_BYTE, 
                    this->cWinImColorData );
      }

    if( this->cValidOverlayData && this->cViewOverlayData ) 
      {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glDrawPixels(this->cWinDataSizeX, this->cWinDataSizeY, GL_RGBA, 
        GL_UNSIGNED_BYTE, this->cWinOverlayData);
      glDisable(GL_BLEND);
      }

    if( this->cViewClickedPoints )
    {
        glColor3f( 0.8, 0.2, 0.6 );
        //glColor3f( 0.8, 0.4, 0.4 );
        glPointSize( 3.0 );
        glBegin(GL_POINTS);
        {
            for ( int i = 0; i < this->numClickedPointsStored(); i++ )
            {
                ClickPoint p;
                this->getClickedPoint( i, p );                
		if (p.clase % 10 == 1) {
		  glColor3f( 1, 1, 1 );
		}
		if (p.clase % 10 == 2) {
		  glColor3f( 0.8, 0.7, 0.2 );
		}
		if (p.clase % 10 == 3) {
		  glColor3f( 0.8, 0.4, 0.4 );
		}
		if (p.clase % 10 == 4) {
		  glColor3f( 1, 0, 0 );
		}
		if (p.clase % 10 == 5) {
		  glColor3f( 0, 1, 0 );
		}		
		if (p.clase % 10 == 6) {
		  glColor3f( 0, 0, 1 );
		}
		if (p.clase % 10 == 7) {
		  glColor3f( 1, 1, 0 );
		}
		if (p.clase % 10 == 8) {
		  glColor3f( 0, 1, 1 );
		}
		if (p.clase % 10 == 9) {
		  glColor3f( 1, 0, 1 );
		}
		if (p.clase % 10 == 0) {
		  glColor3f( 0, 0, 0 );
		}
 
                float pts[3] = { p.x, p.y, p.z };
		
                if ( static_cast<int>( pts[this->cWinOrder[2]] ) ==
                     (int)this->sliceNum() )
                {
                    float x;
                    if(this->cFlipX[this->cWinOrientation])
                    {
                        x = this->cW - (pts[this->cWinOrder[0]] 
                                        - this->cWinMinX) * scale0
                            - originX;
                    }
                    else
                    {
                        x = (pts[this->cWinOrder[0]] - this->cWinMinX) * scale0
                            + originX;
                    }

                    float y;
                    if(this->cFlipY[this->cWinOrientation])
                    {
                        y = this->cH - (pts[this->cWinOrder[1]] 
                                        - this->cWinMinY) * scale1
                            - originY;
                    }
                    else
                    {
                        y = (pts[this->cWinOrder[1]] - this->cWinMinY) * scale1
                             + originY;
                    }
                    glVertex2f( x, y );
                }
            }
        }
        glEnd();
    }

    if( this->cViewAxisLabel ) 
      {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glColor4f(0.2, 0.2, 0.78, (float)0.75);
      gl_font(FL_TIMES_BOLD, 12);
      
      if( !this->cFlipX[this->cWinOrientation] )
        {
        const int y = static_cast<int>(  this->cH/2-gl_height()/2  );
        gl_draw( this->cAxisLabelX[this->cWinOrientation],
          this->cW-(gl_width(this->cAxisLabelX[this->cWinOrientation])+10), 
          static_cast<float>( y ) );
        }
      else
        {
        const int y = static_cast<int>( this->cH/2-gl_height()/2  );
        gl_draw( this->cAxisLabelX[this->cWinOrientation],
          (gl_width(this->cAxisLabelX[this->cWinOrientation])+10),
          static_cast<float>( y ));
        }
      
      if(!this->cFlipY[this->cWinOrientation])
        {
        const int y = static_cast<int>( this->cH-gl_height()-10 ) ;
        gl_draw( this->cAxisLabelY[this->cWinOrientation],
          this->cW/2-(gl_width(this->cAxisLabelY[this->cWinOrientation])/2),
          static_cast<float>(y) );
        }
      else
        {
        const int y = static_cast<int>( gl_height()+10 );
        gl_draw( this->cAxisLabelY[this->cWinOrientation], 
          this->cW/2-(gl_width(this->cAxisLabelY[this->cWinOrientation])/2),
          static_cast<float>(y));
        }
      
      glDisable(GL_BLEND);
      }
    if( this->cViewValue ) 
      {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glColor4f(0.1, 0.64, 0.2, (float)0.75);
      gl_font(FL_TIMES_BOLD, 12);
      char s[80];
      float px, py, pz, val = this->cClickSelectV;
      char * suffix = (char*)"";
      if( this->cViewValuePhysicalUnits )
        {
        px = this->cOrigin[0]+this->cSpacing[0]*this->cClickSelect[0];
        py = this->cOrigin[1]+this->cSpacing[1]*this->cClickSelect[1];
        pz = this->cOrigin[2]+this->cSpacing[2]*this->cClickSelect[2];
        suffix = this->cPhysicalUnitsName;
        }
       else
        {
        px = this->cClickSelect[0];
        py = this->cClickSelect[1];
        pz = this->cClickSelect[2];
        }
      if((ImagePixelType)1.5==1.5)
        {
        sprintf(s, "(%0.1f%s,  %0.1f%s,  %0.1f%s) = %0.3f", 
                px, suffix,
                py, suffix,
                pz, suffix,
                val);
        }
      else
        {
        sprintf(s, "(%0.1f%s,  %0.1f%s,  %0.1f%s) = %d", 
                px, suffix,
                py, suffix,
                pz, suffix,
                (int)val);
        }
      gl_draw( s,
        (int)(this->cW-(gl_width(s)+2)), 2);
      glDisable(GL_BLEND);
      
      }
	if( this->cViewRAStags ) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glColor4f(0.9, 0.2, 0.2, (float)0.75);
      gl_font(FL_TIMES_BOLD, 15);
   	  float v1[3] = {this->cOrientationMatrix[0][0],this->cOrientationMatrix[1][0],this->cOrientationMatrix[2][0]};
	  float v2[3] = {this->cOrientationMatrix[0][1],this->cOrientationMatrix[1][1],this->cOrientationMatrix[2][1]};
	  float v3[3] = {this->cOrientationMatrix[0][2],this->cOrientationMatrix[1][2],this->cOrientationMatrix[2][2]};
	  if (this->cWinOrientation == 0) { // Dir X	    
		 //std::cout << "Dir X. v3=[" << v3[0] << " " << v3[1] << " " << v3[2]<< "]   v2=[" << v2[0] << " " << v2[1] << " " << v2[2] << "]" << std::endl;
	     DrawLateralLetters(v3);
		 DrawVerticalLetters(v2);
	  }
	  if (this->cWinOrientation == 1) { // Dir Y
	     //std::cout << "Dir Y. v1=[" << v1[0] << " " << v1[1] << " " << v1[2]<< "]   v3=[" << v3[0] << " " << v3[1] << " " << v3[2] << "]" << std::endl;
		 DrawLateralLetters(v1);
		 DrawVerticalLetters(v3);
	  }
	  if (this->cWinOrientation == 2) { // Dir Z
		 //std::cout << "Dir Z. v1=[" << v1[0] << " " << v1[1] << " " << v1[2]<< "]   v2=[" << v2[0] << " " << v2[1] << " " << v2[2] << "]" << std::endl;
		 DrawLateralLetters(v1);
		 DrawVerticalLetters(v2);
	  }
	  
	  glDisable(GL_BLEND);
	  
	}
    if( this->cViewDetails )
      {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glColor4f(0.9, 0.4, 0.1, (float)0.75);
      gl_font(FL_TIMES_BOLD, 12);
      char s[80];
      if(this->cWinOrientation == 0)
        sprintf(s, "X - Slice: %3d", this->cWinCenter[0]);
      else if(this->cWinOrientation == 1)
        sprintf(s, "Y - Slice: %3d", this->cWinCenter[1]);
      else
        sprintf(s, "Z - Slice: %3d", this->cWinCenter[2]);
      gl_draw( s, 2, 2+5*(gl_height()+2) );
      sprintf(s, "Dims: %3d x %3d x %3d", 
        (int)this->cDimSize[0], (int)this->cDimSize[1], (int)this->cDimSize[2]);
      gl_draw( s, 2, 2+4*(gl_height()+2) );
      sprintf(s, "Voxel: %0.3f x %0.3f x %0.3f", 
        this->cSpacing[0], this->cSpacing[1], this->cSpacing[2]);
      gl_draw( s, 2, 2+3*(gl_height()+2) );
      sprintf(s, "Int. Range: %0.3f - %0.3f", (float)this->cDataMin, 
              (float)this->cDataMax);
      gl_draw( s, 2, 2+2*(gl_height()+2) );
      sprintf(s, "Int. Window: %0.3f(%s) - %0.3f(%s)", 
        (float)this->cIWMin, IWModeTypeName[this->cIWModeMin], 
        (float)this->cIWMax, IWModeTypeName[this->cIWModeMax]);
      gl_draw( s, 2, 2+1*(gl_height()+2) );
      sprintf(s, "View Mode: %s", ImageModeTypeName[this->cImageMode]);
      gl_draw( s, 2, 2+0*(gl_height()+2) );
      glDisable(GL_BLEND);
      }
    if( this->cViewCrosshairs 
      && static_cast<int>(this->cClickSelect[this->cWinOrder[2]]) == 
         static_cast<int>( this->sliceNum() ) )
      {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glColor4f(0.1, 0.64, 0.2, (float)0.75);
      int x;
      if(this->cFlipX[this->cWinOrientation])
        {
        x = (int)(this->cW - (this->cClickSelect[this->cWinOrder[0]] 
                           - this->cWinMinX) * scale0 - originX);
        }
      else
        {
        x = (int)((this->cClickSelect[this->cWinOrder[0]] 
                   - this->cWinMinX) * scale0 + originX);
        }
      int y;
      if(this->cFlipY[this->cWinOrientation])
        {
        y = (int)(this->cH - (this->cClickSelect[this->cWinOrder[1]] 
                           - this->cWinMinY) * scale1 - originY);
        }
      else
        {
        y = (int)((this->cClickSelect[this->cWinOrder[1]] 
                   - this->cWinMinY) * scale1 + originY);
        }
      glBegin(GL_LINES);
      glVertex2d(0, y);
      glVertex2d(x-2, y);
      glVertex2d(x+2, y);
      glVertex2d(this->cW-1, y);
      glVertex2d(x, 0);
      glVertex2d(x, y-2);
      glVertex2d(x, y+2);
      glVertex2d(x, this->cH-1);
      glEnd();
      glDisable(GL_BLEND);
      }
    }
  }





//
//
template <class ImagePixelType, class OverlayPixelType>
int GLSliceView<ImagePixelType, OverlayPixelType>::handle(int event)
  {
  int x = Fl::event_x();
  int y = Fl::event_y();
  int button;
  
  static int boxX, boxY;
  double zoomBase = this->cW/(this->cDimSize[this->cWinOrder[0]]*(fabs(this->cSpacing[this->cWinOrder[0]])/fabs(this->cSpacing[0])));
	  if(zoomBase >
		 this->cH/(this->cDimSize[this->cWinOrder[1]]*(fabs(this->cSpacing[this->cWinOrder[1]])/fabs(this->cSpacing[0]))))
	  {
		  zoomBase = this->cH/(this->cDimSize[this->cWinOrder[1]]*(fabs(this->cSpacing[this->cWinOrder[1]])/fabs(this->cSpacing[0])));
	  }
	  
  double scale0 = this->cWinZoom * zoomBase * fabs(this->cSpacing[this->cWinOrder[0]])/fabs(this->cSpacing[0]);
  double scale1 = this->cWinZoom * zoomBase * fabs(this->cSpacing[this->cWinOrder[1]])/fabs(this->cSpacing[0]);
	  
  switch(event)
    {      
      case FL_PUSH:
      case FL_DRAG:     
        button = Fl::event_button()-1;               
        if(button <= 0) 
          {
          if(this->cClickMode == CM_BOX || this->cClickMode == CM_SELECT ) 
			{
			  double originX = 0;
			  double originY = 0;
			  if(this->cWinZoom<=1)
				{
				if(this->cW-scale0*this->cDimSize[this->cWinOrder[0]]>0)
				  {
					originX = (int)((this->cW-scale0*this->cDimSize[this->cWinOrder[0]])/2.0);
				  }
					if(this->cH-scale1*this->cDimSize[this->cWinOrder[1]]>0)
				  {
					originY = (int)((this->cH-scale1*this->cDimSize[this->cWinOrder[1]])/2.0);
				  }
				}
				float p[3];
				p[this->cWinOrder[0]] = this->cWinMinX + ( (1-this->cFlipX[this->cWinOrientation])*(x-originX) 
											+ (this->cFlipX[this->cWinOrientation])*(this->cW-x-originX) ) / scale0;
					if(p[this->cWinOrder[0]]<this->cWinMinX) 
						p[this->cWinOrder[0]] = this->cWinMinX;
					if(p[this->cWinOrder[0]]>this->cWinMaxX) 
						p[this->cWinOrder[0]] = this->cWinMaxX;
					p[this->cWinOrder[1]] = this->cWinMinY + (this->cFlipY[this->cWinOrientation]*(y-originY) 
												  + (1-this->cFlipY[this->cWinOrientation])*(this->cH-y-originY)) / scale1;
					if(p[this->cWinOrder[1]]<this->cWinMinY) 
						p[this->cWinOrder[1]] = this->cWinMinY;
					if(p[this->cWinOrder[1]]>this->cWinMaxY) 
						p[this->cWinOrder[1]] = this->cWinMaxY;
					if(this->cImageMode != IMG_MIP)
						p[this->cWinOrder[2]] = this->cWinCenter[this->cWinOrder[2]];
					else
						p[this->cWinOrder[2]] = this->cWinZBuffer[(int)p[this->cWinOrder[0]]
													  - this->cWinMinX 
													  + ((int)p[this->cWinOrder[1]]
														 - this->cWinMinY)
													  * this->cWinDataSizeX];
					if(this->cClickMode == CM_SELECT)
						//std::cout << "clicking in handle GLSliceview\n";
						SetBuble(p[0], p[1], p[2], this->m_radius, this->cActiveClass);
					 
					return 1;
				}
				
            //if(event == FL_PUSH)            
//              {            
//              boxX = x;
//              boxY = y;
//              }
//            else
//              {
//              if(event == FL_DRAG)
//                {	     
//                make_current();
//                fl_overlay_clear();
//                fl_overlay_rect(boxX, boxY, x-boxY, y-boxY);
//                }
//              else
//                {              
//                make_current();
//                fl_overlay_clear();
//                }
//              }
         }
		return 0;	
        break;
	case FL_RELEASE:
    default:
      break;
    }
  
  return SliceView<ImagePixelType>::handle(event);
  }

}; //namespace
#endif

