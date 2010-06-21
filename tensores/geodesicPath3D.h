#ifndef __geodesicPath3D_h
#define __geodesicPath3D_h

#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <itkLevelSet.h>
#include <itkDerivativeImageFilter.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <queue>
#include <vector>

namespace itk
{


template <class TInput, class TOutput>
class ITK_EXPORT geodesicPath3D :
    public ImageToImageFilter< TInput, TOutput > 
{

public:
  /** Standard class typedefs. */
  typedef geodesicPath3D  Self;
  typedef ImageToImageFilter<TInput,TOutput>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(geodesicPath3D, ImageToImageFilter);

  /** Some convenient typedefs. */

  typedef TInput InputType;
  typedef typename InputType::ConstPointer InputPointer;
  typedef typename InputType::RegionType InputRegionType; 
  typedef typename InputType::PixelType InputPixelType; 
  typedef typename InputType::IndexType InputIndexType; 
  typedef typename	InputType::PointType  InputPointType; 
  typedef typename	InputType::SpacingType	InputSpacingType; 
  typedef typename	InputType::SizeType	InputSizeType; 

  typedef TOutput OutputType;
  typedef typename OutputType::Pointer OutputPointer;
  typedef typename OutputType::RegionType OutputRegionType;
  typedef typename OutputType::PixelType OutputPixelType;
  typedef typename OutputType::SizeType OutputSizeType;
  typedef typename OutputType::SpacingType OutputSpacingType;
  typedef typename OutputType::PointType   OutputPointType;
  typedef typename OutputType::IndexType OutputIndexType; 
	  
  itkStaticConstMacro(ImageDimension, unsigned int,TInput::ImageDimension);

  typedef LevelSetTypeDefault<TOutput>  LevelSetType;
  typedef typename LevelSetType::NodeType NodeType;
  typedef typename LevelSetType::NodeContainer NodeContainer; 
  typedef typename LevelSetType::NodeContainer::Pointer NodeContainerPointer;
  

  typedef typename itk::ImageRegionIterator<OutputType> IteratorType;
  typedef typename IteratorType::OffsetType OutOffsetType;
  
  typedef typename itk::DerivativeImageFilter< InputType, InputType > GradFilterType; 
  typedef itk::CovariantVector<float,3>         CovPixelType; 
  typedef itk::Image<CovPixelType, 3 > CovImageType; 
  typedef typename itk::GradientRecursiveGaussianImageFilter< InputType, CovImageType > GradRecursiveGaussianFilterType;

  void SetInput( const TInput * input);
  itkSetMacro(max_num_points,unsigned long);  
  itkSetMacro(Verbose, bool);
  itkGetMacro(max_num_points,unsigned long);
  void SetSeeds(NodeContainer * points);
  void GetNewCoordinates3D( int *newx,int *newy,int* newz,int x,int y,int z, 
                     float gradientx,float gradienty,float gradientz,
					 float angle_fi_error,float angle_theta_error,
					 float *new_fi_error,float *new_theta_error);
  //NodeContainer* GetSeeds();
  typedef typename std::vector < std::vector < OutputIndexType >  > VectorContainerType;
  typedef typename itk::Vector<float,3> VectorFloatType;
  typedef typename std::vector < std::vector < InputPointType >  > VectorContainerFloatType;
	
  typename std::vector < std::vector < OutputIndexType >  >  GetPointsOut() {
     return m_PointsOut;
  }
   
  typename std::vector < std::vector < InputPointType >  >  GetPointsOutFloat() {
     return m_PointsOutFloat;
  }
 
   virtual void SetOutputSize( const OutputSizeType& size )
  { m_OutputRegion = size; }
  virtual OutputSizeType GetOutputSize() const
  { return m_OutputRegion.GetSize(); }
  itkSetMacro( OutputRegion, OutputRegionType );
  itkGetConstReferenceMacro( OutputRegion, OutputRegionType );
  itkSetMacro( OutputSpacing, OutputSpacingType );
  itkGetConstReferenceMacro( OutputSpacing, OutputSpacingType );
  itkSetMacro( OutputOrigin, OutputPointType );
  itkGetConstReferenceMacro( OutputOrigin, OutputPointType );
  itkSetMacro( OverrideOutputInformation, bool );
  itkGetConstReferenceMacro( OverrideOutputInformation, bool );
  itkBooleanMacro( OverrideOutputInformation );

  itkSetMacro( H, float );
  itkSetMacro( Metodo, int );
  itkSetMacro( MetodoPtoIntermedio, int );
  itkSetMacro( Sigma, float );

protected:
  
 geodesicPath3D();
 virtual ~geodesicPath3D() {};
 virtual void GenerateOutputInformation();
 void GenerateData();
 
 void PrintSelf(std::ostream& os, Indent indent) const;
 
private:
 geodesicPath3D(const Self&); //purposely not implemented
 void operator=(const Self&); //purposely not implemented
  
 typedef std::vector<NodeType> HeapContainer;
 typedef std::greater<NodeType> NodeComparer;
 typedef std::priority_queue< NodeType, HeapContainer, NodeComparer > HeapType;
 //typedef std::priority_queue< NodeType> HeapType;
 
 OutputSizeType                                m_OutputSize;
 OutputRegionType                              m_OutputRegion;
 OutputSpacingType                             m_OutputSpacing;
 OutputPointType                               m_OutputOrigin;
 bool                                          m_OverrideOutputInformation;

 HeapType m_Heap;
 unsigned long m_max_num_points;
 NodeContainerPointer m_seeds;
 VectorContainerType m_PointsOut;
 VectorContainerFloatType m_PointsOutFloat;
 bool m_Verbose;
 int *max;
 OutOffsetType v_offset26[26];
 float pi;
 int m_Metodo,m_MetodoPtoIntermedio;
 float m_Sigma,m_H;
  
};

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "geodesicPath3D.txx"
#endif

#endif
