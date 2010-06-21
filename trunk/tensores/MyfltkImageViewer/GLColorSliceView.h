/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: GLColorSliceView.h,v $
 Language:  C++
 Date:      $Date: 2006/05/19 22:31:34 $
 Version:   $Revision: 1.10 $
 
  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
  
   This software is distributed WITHOUT ANY WARRANTY; without even 
   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
   PURPOSE.  See the above copyright notices for more information.
   
=========================================================================*/
#ifndef GLCOLORSLICEVIEW_H
#define GLCOLORSLICEVIEW_H

#include <FL/gl.h>
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Gl_Window.H>

#include "itkColorTable.h"
#include "itkRGBPixel.h"

#include "SliceView.h"
#include "GLSliceView.h"

#include <math.h>

namespace itk {
  
/**
* GLSliceView : Derived from abstract class SliceView and Fl_Gl_Window
* See SliceView.h for details...
  **/
template <class ImagePixelType, class OverlayPixelType>
class GLColorSliceView 
  : public GLSliceView<ImagePixelType, OverlayPixelType> 
  {
  public:

    typedef itk::Image<RGBPixel<ImagePixelType>, 3> ImageType;
    typedef typename ImageType::Pointer ImagePointer;
    typedef typename ImageType::ConstPointer ImageConstPointer;
    typedef GLSliceView<ImagePixelType, OverlayPixelType> Superclass;
    typedef typename Superclass::RegionType RegionType;
    typedef typename Superclass::SizeType SizeType;
    typedef typename Superclass::IndexType IndexType;

    /*! FLTK required constructor - must use imData() to complete 
     *  definition */
    GLColorSliceView(int x, int y, int w, int h, const char *l);
    
    /*! Specify the 3D image to view slice by slice */
    virtual void SetInputImage(ImageType * newImData);
    virtual const ImagePointer & GetInputImage(void);

    virtual void clickSelect(float newX, float newY, float newZ);
    virtual void size(int w, int h);
    
    virtual void resize(int x, int y, int w, int h);

    virtual void update();

    virtual void draw();
    
    virtual int handle(int event);

  protected:
     void clickMoveRGB(float newX, float newY, float newZ);

    float         cClickSelectR;
    float         cClickSelectG;
    float         cClickSelectB;

    ImagePointer  cImData;
  };
  
template <class ImagePixelType, class OverlayPixelType>
GLColorSliceView<ImagePixelType, OverlayPixelType>::
GLColorSliceView(int x, int y, int w, int h, const char *l):
  GLSliceView<ImagePixelType, OverlayPixelType>(x, y, w, h, l)
  {
  cClickSelectR = 0;
  cClickSelectG = 0;
  cClickSelectB = 0;
  }
  
template <class ImagePixelType, class OverlayPixelType>
void GLColorSliceView<ImagePixelType, OverlayPixelType>::
clickMoveRGB(float newX, float newY, float newZ)
  {   
    this->cClickSelect[0] = newX;
  if(this->cClickSelect[0]<0)
    this->cClickSelect[0] = 0;
  if(this->cClickSelect[0] >= this->cDimSize[0])
    this->cClickSelect[0] = this->cDimSize[0]-1;
  
  this->cClickSelect[1] = newY;
  if(this->cClickSelect[1]<0)
    this->cClickSelect[1] = 0;
  if(this->cClickSelect[1] >= this->cDimSize[1])
    this->cClickSelect[1] = this->cDimSize[1]-1;
  
  this->cClickSelect[2] = newZ;
  if(this->cClickSelect[2]<0)
    this->cClickSelect[2] = 0;
  if(this->cClickSelect[2] >= this->cDimSize[2])
    this->cClickSelect[2] = this->cDimSize[2]-1;

  //if(this->cClickMoveCallBack != NULL)
//    this->cClickMoveCallBack(cClickSelect[0], this->cClickSelect[1], 
//    this->cClickSelect[2]);
//  if(this->cClickMoveArgCallBack != NULL)
//    this->cClickMoveArgCallBack(cClickSelect[0], this->cClickSelect[1], 
//    this->cClickSelect[2],cClickMoveArg);

  typename ImageType::IndexType ind;
  
  ind[0] = (unsigned long)this->cClickSelect[0];
  ind[1] = (unsigned long)this->cClickSelect[1];
  ind[2] = (unsigned long)this->cClickSelect[2];
  
 
  cClickSelectR = this->cImData->GetPixel(ind)[0];
  cClickSelectG = this->cImData->GetPixel(ind)[1];
  cClickSelectB = this->cImData->GetPixel(ind)[2];
  this->update();

}
  
//
//
template <class ImagePixelType, class OverlayPixelType>
const typename Image<RGBPixel<ImagePixelType>,3>::Pointer &
GLColorSliceView<ImagePixelType, OverlayPixelType>
::GetInputImage(void)
  {
  return this->cImData;
  }


template <class ImagePixelType, class OverlayPixelType>
void GLColorSliceView<ImagePixelType, OverlayPixelType>::
SetInputImage(ImageType * newImData)
  {
  RegionType region = newImData->GetLargestPossibleRegion();
  if( region.GetNumberOfPixels() == 0 ) 
    {
    return;
    }

  SizeType   size   = region.GetSize();


  // If the overlay has been set and the size is different from the new image,
  // it is removed.
  if( this->cValidOverlayData)
  {  
    SizeType  overlay_size   = this->cOverlayData->GetLargestPossibleRegion().GetSize();
      
    if((overlay_size[0] != size[0])
       ||  (overlay_size[1] != size[1])
       ||  (overlay_size[2] != size[2])
      )
      {
       if(this->cWinOverlayData != NULL)
        {
        delete [] this->cWinOverlayData;
        }
       this->cWinOverlayData       = NULL;
       this->cValidOverlayData     = false;
      }       
  }
    
  this->cImData = newImData;
  this->cDimSize[0] = size[0];
  this->cDimSize[1] = size[1];
  this->cDimSize[2] = size[2];
  this->cSpacing[0] = this->cImData->GetSpacing()[0];
  this->cSpacing[1] = this->cImData->GetSpacing()[1];
  this->cSpacing[2] = this->cImData->GetSpacing()[2];

  //calculating cDataMax and cDataMin
  IndexType ind;
  ind[0] = 0; 
  ind[1] = 0; 
  ind[2] = 0;
  
  this->cDataMax = this->cImData->GetPixel(ind)[0];
  this->cDataMin = this->cDataMax;
  RGBPixel<ImagePixelType> tfv;
  double tf;
    
    
  for( unsigned int i=0; i<this->cDimSize[0]; i++ )
    {
    ind[0] = i;
    for(unsigned int j=0; j<this->cDimSize[1]; j++ )
      {
      ind[1] = j;
      for( unsigned int k=0; k<this->cDimSize[2]; k++ )
        {
        ind[2] = k;
        tfv = this->cImData->GetPixel(ind);
        for( unsigned int l=0; l<3; l++)
          {
          tf = (double)(tfv[l]);
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
    }
  
  this->cIWMin      = this->cDataMin;
  this->cIWMax      = this->cDataMax;
  this->cIWModeMin  = IW_MIN;
  this->cIWModeMax  = IW_MAX;
    
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
  
  this->cWinDataSizeX = this->cDimSize[0];
  this->cWinDataSizeY = this->cDimSize[1];
  
  if(this->cWinImData != NULL)
    {
    delete [] this->cWinImData;
    }
  
  this->cWinImData = new unsigned char[ this->cWinDataSizeX * this->cWinDataSizeY * 3 ];
  
  if(this->cWinZBuffer != NULL) 
    {
    delete [] this->cWinZBuffer;
    }
  
  this->cWinZBuffer = new unsigned short[ this->cWinDataSizeX * this->cWinDataSizeY ];
  
  this->cViewImData  = true;
  this->cValidImData = true;
    
  update();
  }

//
//
template <class ImagePixelType, class OverlayPixelType>
void GLColorSliceView<ImagePixelType, OverlayPixelType>::
clickSelect(float newX, float newY, float newZ)
  {    
  this->cClickSelect[0] = newX;
  if(this->cClickSelect[0]<0) 
    {
    this->cClickSelect[0] = 0;
    }
  if(this->cClickSelect[0] >= this->cDimSize[0])
    {
    this->cClickSelect[0] = this->cDimSize[0]-1;
    }

  this->cClickSelect[1] = newY;
  if(this->cClickSelect[1]<0)
    {
    this->cClickSelect[1] = 0;
    }
  if(this->cClickSelect[1] >= this->cDimSize[1])
    {
    this->cClickSelect[1] = this->cDimSize[1]-1;
    }

  this->cClickSelect[2] = newZ;

  if(this->cClickSelect[2]<0)
    {
    this->cClickSelect[2] = 0;
    }
  if(this->cClickSelect[2] >= this->cDimSize[2])
    {
    this->cClickSelect[2] = this->cDimSize[2]-1;
    }

  typename ImageType::IndexType ind;

  ind[0] = (unsigned long)this->cClickSelect[0];
  ind[1] = (unsigned long)this->cClickSelect[1];
  ind[2] = (unsigned long)this->cClickSelect[2];
  this->cClickSelectV = this->cImData->GetPixel(ind)[0];
  cClickSelectR = cImData->GetPixel(ind)[0];
  cClickSelectG = cImData->GetPixel(ind)[1];
  cClickSelectB = cImData->GetPixel(ind)[2];
      
  /*if length of list is equal to max, remove the earliest point stored */
  if((this->maxClickPoints>0)&&(this->cClickedPoints.size() == this->maxClickPoints))
    {
    this->cClickedPoints.pop_back();
    }
  
  ClickPoint point( this->cClickSelect[0], 
                    this->cClickSelect[1], 
                    this->cClickSelect[2], 
                    this->cClickSelectV, 0 );

  this->cClickedPoints.push_front( point );

  if(this->cClickSelectCallBack != NULL)
    {
      this->cClickSelectCallBack(this->cClickSelect[0], this->cClickSelect[1], 
                           this->cClickSelect[2], this->cClickSelectV, 0);
    }
  if(this->cClickSelectArgCallBack != NULL)
    {
    this->cClickSelectArgCallBack(this->cClickSelect[0], this->cClickSelect[1], 
                              this->cClickSelect[2], this->cClickSelectV, 0, 
                              this->cClickSelectArg);
    }

  if(this->cViewValue || this->cViewCrosshairs)
    {
    this->redraw();
    }
  }


template <class ImagePixelType, class OverlayPixelType>
void GLColorSliceView<ImagePixelType, OverlayPixelType>::
size(int w, int h)
  {
  GLSliceView<ImagePixelType, OverlayPixelType>::size(w, h);
  }
    
template <class ImagePixelType, class OverlayPixelType>
void GLColorSliceView<ImagePixelType, OverlayPixelType>::
resize(int x, int y, int w, int h)
  {
  GLSliceView<ImagePixelType, OverlayPixelType>::resize(x, y, w, h);
  }

//
template <class ImagePixelType, class OverlayPixelType>
void GLColorSliceView<ImagePixelType, OverlayPixelType>::
update()
  {
  if( !this->cValidImData ) 
    {
    return;
    }
	
  // Copiado del GLSliceView 
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
  // Copiado del GLSliceView 
  /////////////////////////  EJE  X ////////////////////////////////////
  if(this->cWinZoom>1) {
    int winWidth = (int)( this->cDimSize[ this->cWinOrder[0] ] / this->cWinZoom );
    //this->cWinSizeX = ( (int) winWidth);
    this->cWinSizeX = (int)( this->cW / scale0 );
    this->cWinMinX = (int)( (int)this->cWinCenter[ this->cWinOrder[0] ] 
                            - this->cWinSizeX/2 );
    this->cWinMaxX = (int)( (int)this->cWinCenter[ this->cWinOrder[0] ] 
                            + this->cWinSizeX/2 );
  } else {
    this->cWinSizeX = (int)(this->cDimSize[ this->cWinOrder[0] ]);
    this->cWinMinX = 0;
    this->cWinMaxX = (int)( (int)(this->cDimSize[ this->cWinOrder[0] ]) - 1 );
    this->cWinCenter[this->cWinOrder[0]] = 
                     (int)( this->cDimSize[ this->cWinOrder[0] ] / 2);
  }
  //int ti = (int)( (int)this->cWinCenter[ this->cWinOrder[0] ] - winWidth/2);
  int ti = this->cWinMinX;
  if( ti <= - (int) this->cDimSize[ this->cWinOrder[0] ] ) 
    {
    ti = -(int)this->cDimSize[ this->cWinOrder[0] ] + 1;
    }
  else if( ti >= (int)this->cDimSize[ this->cWinOrder[0] ]) 
    {
    ti = this->cDimSize[ this->cWinOrder[0] ] - 1;
    }
  this->cWinMinX = ti;
  this->cWinMaxX = this->cDimSize[ this->cWinOrder[0] ] - 1; // here
  if( this->cWinMaxX >= static_cast<int>( this->cDimSize[ this->cWinOrder[0] ] ) )
    {
    this->cWinMaxX = this->cDimSize[ this->cWinOrder[0] ] - 1;
    }
  /////////////////////////  EJE  Y ////////////////////////////////////
  if(this->cWinZoom>1) {
    this->cWinSizeY = (int)( this->cH / scale1 );
    this->cWinMinY = (int)( (int)(this->cWinCenter[ this->cWinOrder[1] ]) 
                             - this->cWinSizeY/2 );
    this->cWinMaxY = (int)( (int)(this->cWinCenter[ this->cWinOrder[1] ]) 
                             + this->cWinSizeY/2 );
  } else {
    this->cWinSizeY = (int)(this->cDimSize[ this->cWinOrder[1] ]);
    this->cWinMinY = 0;
    this->cWinMaxY = (int)( (int)(this->cDimSize[ this->cWinOrder[1] ]) - 1 );
    this->cWinCenter[this->cWinOrder[1]] = 
                     (int)( this->cDimSize[ this->cWinOrder[1] ] / 2);
  }
  //ti = static_cast<int>( static_cast<int>(this->cWinCenter[ this->cWinOrder[1] ]) 
  //                       - winWidth/2);
  ti = this->cWinMinY;
  if( ti <= - static_cast<int>( this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    ti = -(int)this->cDimSize[ this->cWinOrder[1] ] + 1;
    }
  else if( ti >= static_cast<int>(this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    ti = this->cDimSize[ this->cWinOrder[1] ] - 1;
    } 
  this->cWinMinY = ti;
  this->cWinMaxY = this->cDimSize[ this->cWinOrder[1] ] - 1;
  if( this->cWinMaxY >= static_cast<int>( this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    this->cWinMaxY = this->cDimSize[ this->cWinOrder[1] ] - 1;
    }
	
//  std::cout << "cW " << this->cW << " cH " << this->cH << std::endl;
//  std::cout << " cWinSizeX " << this->cWinSizeX << " cWinSizeY " << this->cWinSizeY << std::endl;
//  std::cout << " cWinMinX " << this->cWinMinX << " cWinMaxX " << this->cWinMaxX << std::endl;
//  std::cout << " cWinMinY " << this->cWinMinY << " cWinMaxY " << this->cWinMaxY << std::endl;
//  std::cout << " cWinZoom " << this->cWinZoom << std::endl;
//  std::cout << " scale0 " << scale0 << " scale1 " << scale1 << std::endl;
//  std::cout << "cWinCenter 0 " << this->cWinCenter[this->cWinOrder[0]] << " cWinCenter 1 " << this->cWinCenter[this->cWinOrder[1]] << std::endl;

  memset( this->cWinImData, 0, this->cWinDataSizeX*this->cWinDataSizeY*3 );
  if( this->cValidOverlayData ) 
    {
    memset(this->cWinOverlayData, 0, this->cWinDataSizeX*this->cWinDataSizeY*4);
    }
  
  IndexType ind;
  
  int l, m;
  
  itk::RGBPixel<ImagePixelType> tfv;
  float tf[3];
  
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

      tfv = this->cImData->GetPixel(ind);
      switch( this->cImageMode ) 
        {
        default:
        case IMG_VAL:
          tf[0] = (float)(((float)tfv[0]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          tf[1] = (float)(((float)tfv[1]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          tf[2] = (float)(((float)tfv[2]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          break;
        case IMG_INV:
          tf[0] = (float)((this->cIWMax-(float)tfv[0])/(this->cIWMax-this->cIWMin)*255);
          tf[1] = (float)((this->cIWMax-(float)tfv[1])/(this->cIWMax-this->cIWMin)*255);
          tf[2] = (float)((this->cIWMax-(float)tfv[2])/(this->cIWMax-this->cIWMin)*255);
          break;
        case IMG_LOG:
          tf[0] = (float)(log((float)tfv[0]-this->cIWMin+0.00000001)
                          /log(this->cIWMax-this->cIWMin+0.00000001)*255);
          tf[1] = (float)(log((float)tfv[1]-this->cIWMin+0.00000001)
                          /log(this->cIWMax-this->cIWMin+0.00000001)*255);
          tf[2] = (float)(log((float)tfv[2]-this->cIWMin+0.00000001)
                          /log(this->cIWMax-this->cIWMin+0.00000001)*255);
          break;
        case IMG_DX:
          if(ind[0]>0) 
            {
            tf[0] = (float)(((float)tfv[0]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[1] = (float)(((float)tfv[1]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[2] = (float)(((float)tfv[2]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            ind[0]--;
            tfv = this->cImData->GetPixel(ind);
            tf[0] -= (float)(((float)tfv[0]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[1] -= (float)(((float)tfv[1]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[2] -= (float)(((float)tfv[2]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            ind[0]++;
            tf[0] += 128;
            tf[1] += 128;
            tf[2] += 128;
            } 
          else
            {
            tf[0] = 128;
            tf[1] = 128;
            tf[2] = 128;
            }
          break;
        case IMG_DY:
          if(ind[1]>0) 
            {
            tf[0] = (float)(((float)tfv[0]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[1] = (float)(((float)tfv[1]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[2] = (float)(((float)tfv[2]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            ind[1]--;
            tfv = this->cImData->GetPixel(ind);
            tf[0] -= (float)(((float)tfv[0]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[1] -= (float)(((float)tfv[1]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[2] -= (float)(((float)tfv[2]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            ind[1]++;
            tf[0] += 128;
            tf[1] += 128;
            tf[2] += 128;
            }
          else
            {
            tf[0] = 128;
            tf[1] = 128;
            tf[2] = 128;
            }
          break;
        case IMG_DZ:
          if(ind[2]>0) 
            {
            tf[0] = (float)(((float)tfv[0]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[1] = (float)(((float)tfv[1]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[2] = (float)(((float)tfv[2]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            ind[2]--;
            tfv = this->cImData->GetPixel(ind);
            tf[0] -= (float)(((float)tfv[0]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[1] -= (float)(((float)tfv[1]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            tf[2] -= (float)(((float)tfv[2]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
            ind[2]++;
            tf[0] += 128;
            tf[1] += 128;
            tf[2] += 128;
            }
          else
            {
            tf[0] = 128;
            tf[1] = 128;
            tf[2] = 128;
            }
          break;
        case IMG_BLEND:
          {
          const int tempval = (int)this->cWinCenter[this->cWinOrder[2]]-1;
          int tmpI = ind[this->cWinOrder[2]];
          ind[this->cWinOrder[2]] = (tempval < 0 ) ? 0 : tempval;
          tfv = this->cImData->GetPixel(ind);
          tf[0] = (float)(tfv[0]);
          tf[1] = (float)(tfv[1]);
          tf[2] = (float)(tfv[2]);
          
          ind[this->cWinOrder[2]] = this->cWinCenter[this->cWinOrder[2]];
          tfv = this->cImData->GetPixel(ind);
          tf[0] += (float)(tfv[0])*2;
          tf[1] += (float)(tfv[1])*2;
          tf[2] += (float)(tfv[2])*2;
          
          const int tempval1 = (int)this->cDimSize[this->cWinOrder[2]]-1;
          const int tempval2 = (int)this->cWinCenter[this->cWinOrder[2]]+1;
          ind[this->cWinOrder[2]] = (tempval1 < tempval2 ) ? tempval1 : tempval2;
          tfv = this->cImData->GetPixel(ind);
          tf[0] += (float)(tfv[0]);
          tf[1] += (float)(tfv[1]);
          tf[2] += (float)(tfv[2]);
          
          tf[0] = (float)((tf[0]/4-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          tf[1] = (float)((tf[1]/4-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          tf[2] = (float)((tf[2]/4-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          ind[this->cWinOrder[2]] = tmpI;
          break;
          }
        case IMG_MIP:
          tf[0] = this->cIWMin;
          tf[1] = this->cIWMin;
          tf[2] = this->cIWMin;
          m = (j-this->cWinMinX) + (k-this->cWinMinY)*this->cWinDataSizeX;
          this->cWinZBuffer[m] = 0;
          int tmpI = ind[this->cWinOrder[2]];
          float tfp = 0;
          float tft = tf[0]+tf[1]+tf[2];
          for(l=0; l<(int)this->cDimSize[this->cWinOrder[2]]; l++) 
            {
            ind[this->cWinOrder[2]] = l;        
            tfv = this->cImData->GetPixel(ind);
            tfp = (float)tfv[0] 
                  + tfv[1] 
                  + tfv[2];
            if(tfp > tft) 
              {
              tf[0] = (float)(tfv[0]);
              tf[1] = (float)(tfv[1]);
              tf[2] = (float)(tfv[2]);
              this->cWinZBuffer[m] = (unsigned short)l;
              }
            }
          tf[0] = (float)((tf[0]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          tf[1] = (float)((tf[1]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          tf[2] = (float)((tf[2]-this->cIWMin)/(this->cIWMax-this->cIWMin)*255);
          ind[this->cWinOrder[2]] = tmpI;
          break;
        }
      
      int c;
      for(c=0; c<3; c++)
        {
        if( tf[c] > 255 )
          {
          switch(this->cIWModeMax) 
            {
            case IW_MIN:
              tf[c] = 0;
              break;
            default:
            case IW_MAX:
              tf[c] = 255;
              break;
            case IW_FLIP:
              tf[c] = 512-tf[c];
              if(tf[c]<0) 
                {
                tf[c] = 0;
                }
              break;
            }
          }
        else 
          {
          if( tf[c] < 0 )
            {
            switch(this->cIWModeMin) 
              {
              default:
              case IW_MIN:
                tf[c] = 0;
                break;
              case IW_MAX:
                tf[c] = 255;
                break;
              case IW_FLIP:
                tf[c] = -tf[c];
                if(tf[c]>255)
                  {
                  tf[c] = 255;
                  }
                break;
              }
            }
          }
        }
        
      l = (j-this->cWinMinX)*3 + (k-this->cWinMinY)*this->cWinDataSizeX*3;
      this->cWinImData[l+0] = (unsigned char)tf[0];
      this->cWinImData[l+1] = (unsigned char)tf[1];
      this->cWinImData[l+2] = (unsigned char)tf[2];
      
      if( this->cValidOverlayData ) 
        {
        l = (j-this->cWinMinX)*4 + (k-this->cWinMinY)*this->cWinDataSizeX*4;
        if(this->cImageMode == IMG_MIP)
          {
          ind[this->cWinOrder[2]] = this->cWinZBuffer[(j-this->cWinMinX) + 
                              (k-this->cWinMinY)*this->cWinDataSizeX];
          }
        
        if( sizeof( OverlayPixelType ) == 1 )
          {
          m = (int)*((unsigned char *)&(this->cOverlayData->GetPixel(ind)));
          if( m > 0 ) {
            m = m - 1;
            this->cWinOverlayData[l+0] = 
              (unsigned char)(this->cColorTable->GetColorComponent(m, 'r')*255);
            this->cWinOverlayData[l+1] = 
              (unsigned char)(this->cColorTable->GetColorComponent(m, 'g')*255);
            this->cWinOverlayData[l+2] = 
              (unsigned char)(this->cColorTable->GetColorComponent(m, 'b')*255);
            this->cWinOverlayData[l+3] = 
              (unsigned char)(this->cOverlayOpacity*255);
            }
          }
        else 
          {
          if(((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[0]
            + ((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[1]
            + ((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[2] > 0)
            {
            if( sizeof( OverlayPixelType ) == 3 )
              {
              this->cWinOverlayData[l+0] = 
                ((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[0];
              this->cWinOverlayData[l+1] = 
                ((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[1];
              this->cWinOverlayData[l+2] = 
                ((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[2];
              this->cWinOverlayData[l+3] = 
                (unsigned char)(this->cOverlayOpacity*255);
              }
            else 
              {
              if( sizeof( OverlayPixelType ) == 4 ) 
                {
                this->cWinOverlayData[l+0] = 
                  ((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[0];
                this->cWinOverlayData[l+1] = 
                  ((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[1];
                this->cWinOverlayData[l+2] = 
                  ((unsigned char *)&(this->cOverlayData->GetPixel(ind)))[2];
                this->cWinOverlayData[l+3] = 
                  (unsigned char)(((unsigned char *)
                  &(this->cOverlayData->GetPixel(ind)))[3]*this->cOverlayOpacity);
                }
              }
            }
          }
        }
      }
    }
  
  this->redraw();
  
  }

//
//
template <class ImagePixelType, class OverlayPixelType>
void GLColorSliceView<ImagePixelType, OverlayPixelType>::
draw()
  {
  if( !this->valid() )
    {
    glClearColor((float)0.0, (float)0.0, (float)0.0, (float)0.0);          
    glShadeModel(GL_FLAT);
    glClear(GL_COLOR_BUFFER_BIT);    //this clears and paints to black
    
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);  //if you don't include this
    //image size differences distort
    //glPixelStorei(GL_PACK_ALIGNMENT, 1);
    }
  else
    {
    glClearColor((float)0.0, (float)0.0, (float)0.0, (float)0.0);          
    glShadeModel(GL_FLAT);
    glClear(GL_COLOR_BUFFER_BIT);    //this clears and paints to black
    
    glMatrixMode(GL_MODELVIEW);    //clear previous 3D draw params
    glLoadIdentity();
    
    glMatrixMode(GL_PROJECTION);
    this->ortho();
    
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

    //float scale0 = this->cW/(float)this->cDimSize[0] * this->cWinZoom
    //  * fabs(this->cSpacing[this->cWinOrder[0]])/fabs(this->cSpacing[0]);
    //float scale1 = this->cW/(float)this->cDimSize[0] * this->cWinZoom
    //  * fabs(this->cSpacing[this->cWinOrder[1]])/fabs(this->cSpacing[0]);
    
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

    //glRasterPos2i((this->cFlipX[this->cWinOrientation])?this->cW:0,
    //  (this->cFlipY[this->cWinOrientation])?this->cH:0);  

    glRasterPos2i((this->cFlipX[this->cWinOrientation])?this->cW-originX:originX,
      (this->cFlipY[this->cWinOrientation])?this->cH-originY:originY);  
    glPixelZoom((this->cFlipX[this->cWinOrientation])?-scale0:scale0,
      (this->cFlipY[this->cWinOrientation])?-scale1:scale1);
    
    if( this->cValidImData && this->cViewImData )
      {
      glDrawPixels( this->cWinDataSizeX, this->cWinDataSizeY, 
                    GL_RGB, GL_UNSIGNED_BYTE, 
                    this->cWinImData );
      }
    
    if( this->cValidOverlayData && this->cViewOverlayData ) 
      {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glDrawPixels(this->cWinDataSizeX, this->cWinDataSizeY, GL_RGBA, 
        GL_UNSIGNED_BYTE, this->cWinOverlayData);
      glDisable(GL_BLEND);
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
      //if((ImagePixelType)1.1==1.1)
      if ( 1 )  
		{
        sprintf(s, "(%0.1f,  %0.1f,  %0.1f) = %0.3f,%03f,%0.3f", 
          this->cClickSelect[0],
          this->cClickSelect[1], 
          this->cClickSelect[2], 
          (float)cClickSelectR,
          (float)cClickSelectG,
          (float)cClickSelectB);
        }
      else
        {
        sprintf(s, "(%0.1f,  %0.1f,  %0.1f) = %d,%d,%d", 
          this->cClickSelect[0],
          this->cClickSelect[1], 
          this->cClickSelect[2], 
          (int)cClickSelectR,
          (int)cClickSelectG,
          (int)cClickSelectB);
        }
      gl_draw( s,
        (int)(this->cW-(gl_width(s)+2)), 2);
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
      //int x;
//      if(this->cFlipX[this->cWinOrientation])
//        {
//        x = (int)(this->cW - (this->cClickSelect[this->cWinOrder[0]] - this->cWinMinX) * scale0);
//        }
//      else
//        {
//        x = (int)((this->cClickSelect[this->cWinOrder[0]] - this->cWinMinX) * scale0);
//        }
//      int y;
//      if(this->cFlipY[this->cWinOrientation])
//        {
//        y = (int)(this->cH - (this->cClickSelect[this->cWinOrder[1]] - this->cWinMinY) * scale1);
//        }
//      else
//        {
//        y = (int)((this->cClickSelect[this->cWinOrder[1]] - this->cWinMinY) * scale1);
//        }
		///////////////////////////////////////////////////////
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

template <class ImagePixelType, class OverlayPixelType>
int GLColorSliceView<ImagePixelType, OverlayPixelType>::
handle(int event)
  {
  
  int x = Fl::event_x();
  int y = Fl::event_y();
  int button;
  
  static int boxX, boxY;
  double zoomBase = this->cW/(this->cDimSize[this->cWinOrder[0]]*(fabs(this->cSpacing[this->cWinOrder[0]])/fabs(this->cSpacing[0])));
  if(zoomBase > this->cH / (this->cDimSize[this->cWinOrder[1]]
                            * (fabs(this->cSpacing[this->cWinOrder[1]])
                              / fabs(this->cSpacing[0]))))
    {
    zoomBase = this->cH / (this->cDimSize[this->cWinOrder[1]]
                           * (fabs(this->cSpacing[this->cWinOrder[1]])
                              / fabs(this->cSpacing[0])));
    }

  double scale0 = this->cWinZoom * zoomBase * fabs(this->cSpacing[this->cWinOrder[0]])/fabs(this->cSpacing[0]);
  double scale1 = this->cWinZoom * zoomBase * fabs(this->cSpacing[this->cWinOrder[1]])/fabs(this->cSpacing[0]);
  
  switch(event) {
  case FL_MOVE:
      button = Fl::event_button()-1;
	
      if(button <= 0)
        {
        if(this->cClickMode == CM_SELECT || this->cClickMode == CM_RGB) 
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
            + (this->cFlipX[this->cWinOrientation])*(this->cW-x-originX) ) 
            / scale0;
          if(p[this->cWinOrder[0]]<this->cWinMinX) 
            p[this->cWinOrder[0]] = this->cWinMinX;
          if(p[this->cWinOrder[0]]>this->cWinMaxX) 
            p[this->cWinOrder[0]] = this->cWinMaxX;
          p[this->cWinOrder[1]] = this->cWinMinY + (this->cFlipY[this->cWinOrientation]*(y-originY) 
            + (1-this->cFlipY[this->cWinOrientation])*(this->cH-y-originY)) 
            / scale1;
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
          if(this->cClickMode == CM_SELECT || this->cClickMode == CM_RGB)  {
            this->clickMoveRGB(p[0], p[1], p[2]);
	      } 
	    }
	  }
      return 0;
      break;   
   default:
      return GLSliceView<ImagePixelType, OverlayPixelType>::handle(event);
  }
}

}; //namespace
#endif

