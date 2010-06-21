/*=========================================
 Copyright (c) Ruben Cardenes Almeida. August 2006
 This work is licensed under a 
 Creative Commons Atribution 2.5 License 
 http://creativecommons.org/licenses/by/2.5
 ==========================================*/

#ifndef __classifyKNNCore_h
#define __classifyKNNCore_h

#include "itkImageToImageFilter.h"
#include <itkImage.h>

namespace itk {

template <class TCh1, class TOutput, unsigned int m_ImageDimension>
class ITK_EXPORT classifyKNNCore :
   public ImageToImageFilter< TCh1, TOutput> 
{
public:
  typedef TOutput OutputType;
  typedef typename OutputType::Pointer OutputPointer;
  typedef typename OutputType::RegionType OutputRegionType; 
  typedef typename OutputType::PixelType OutputPixelType; 
  typedef typename OutputType::SizeType OutputSizeType;
  typedef TCh1 Ch1Type;
  typedef typename Ch1Type::Pointer Ch1Pointer;
  typedef typename Ch1Type::RegionType Ch1RegionType; 
  typedef typename Ch1Type::PixelType Ch1PixelType; 
  typedef typename Ch1Type::SizeType Ch1SizeType;
  typedef typename OutputType::IndexType IndexType;
  typedef itk::Image<std::vector<int>, 1> DiagramImageType;
  typedef typename DiagramImageType::Pointer DiagramPointer;
  typedef typename DiagramImageType::RegionType DiagramRegionType; 
  typedef typename DiagramImageType::PixelType DiagramPixelType; 
  typedef typename DiagramImageType::SizeType DiagramSizeType;
  typedef typename DiagramImageType::IndexType DiagramIndexType;

  typedef classifyKNNCore  Self;
  typedef ImageToImageFilter<TCh1,TOutput>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

 /** Run-time type information (and related methods). */
  itkTypeMacro(classifyKNNCore, ImageToImageFilter);

  void SetCh1( const TCh1 * channel1);
  void SetDiagram( DiagramImageType * diagram);

  struct NodeType {
    DiagramIndexType index_intensity;
    IndexType index_space;
    int clase;
  };

  typedef  std::vector<NodeType> ListType;

  void SetPrototipos (const ListType _arg) 
   { 
    //Luego lo volvere a volcar en Lista 0, pero es un mal menor.
    //Later, we will put this in List 0
     m_Prototipos.insert(m_Prototipos.begin(),_arg.begin(),_arg.end());
       this->Modified(); 
   } 
 
 itkSetMacro(K, unsigned short);

protected:

  classifyKNNCore();
  virtual ~classifyKNNCore() {};

  void GenerateOutputInformation();
  void GenerateData();
  // void PrintSelf(std::ostream& os, Indent indent) const;
 
  private:
  classifyKNNCore(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //unsigned short m_MaxClass;
  ListType m_Prototipos;
  DiagramPointer voro;
  unsigned short m_K;
};


}//end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
#include "classifyKNNCore.txx"
#endif

#endif

