/*
 *  itkTractography.h
 *
 */
 
#ifndef __itkTractography_h
#define __itkTractography_h

#include "itkImageToImageFilter.h"

namespace itk
{

template <class TInput, class TOutput>
class ITK_EXPORT Tractography: public ImageToImageFilter<TInput, TOutput>
{
public:
  /** Standard class typedefs. */
  typedef Tractography  Self;
  typedef ImageToImageFilter<TInput,TOutput>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Tractography, ImageToImageFilter);

 /** Some convenient typedefs. */
   typedef TInput													InputImageType;
	typedef TOutput													OutputImageType;
    typedef typename InputImageType::Pointer                        InputImagePointer;
	typedef typename InputImageType::ConstPointer					InputImageConstPointer;
	typedef typename OutputImageType::Pointer                       OutputImagePointer;
	typedef typename InputImageType::RegionType						InputImageRegionType;
	
    typedef typename		InputImageType::Pointer			InputPointer;
	typedef typename	InputImageType::IndexType		InputIndexType; 
	typedef typename	InputImageType::PointType		InputPointType; 
	typedef typename	InputImageType::SpacingType		InputSpacingType; 
	
	typedef typename    itk::Image<unsigned char, 3>	CharImageType;
	typedef typename    CharImageType::Pointer			CharImageTypePointer;
	typedef typename	CharImageType::IndexType		CharIndexType; 
	typedef typename	CharImageType::PointType		CharPointType; 
	typedef typename	CharImageType::SpacingType		CharSpacingType; 

	typedef typename	itk::Vector<InputIndexType>		IndexVectorType;
	
/*    typedef typename	OutputImageType::Pointer		OutputPointer;*/
	typedef typename	OutputImageType::RegionType		OutputImageRegionType;
/*	typedef typename	OutputImageType::PixelType		OutputPixelType;*/
	typedef typename	OutputImageType::SizeType		OutputSizeType;
/*	typedef typename	OutputImageType::SpacingType	OutputSpacingType;
	typedef typename	OutputImageType::PointType		OutputPointType;
	typedef typename	OutputImageType::IndexType		OutputIndexType; */
	  
  typedef typename std::vector<InputPointType>						StreamlineType;
  typedef typename std::vector<StreamlineType>						StreamlineListType;

  itkStaticConstMacro(ImageDimension, unsigned int,TInput::ImageDimension);
  
  itkSetMacro(CurvatureThreshold, float);
  itkSetMacro(StepLength, float);
  itkSetMacro(FaThreshold, float);
  itkSetMacro(MinimumFiberLength, unsigned int);
  itkSetMacro(Connection, bool);

  itkSetMacro(FibersPerVoxel, unsigned int);
  itkGetMacro(FibersPerVoxel, unsigned int);
  void GetStreamlineList(StreamlineListType& list){
	list=m_StreamlineList;
	}
	
  itkSetMacro(SeedsImage, CharImageTypePointer);
	
 Tractography();
 virtual ~Tractography() {};
 virtual void GenerateOutputInformation();
  //void GenerateData();
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId);
  void RungeKuttaIntegrationWithConnection( );
  void RungeKuttaIntegration(); 
  void RungeKuttaIntegration(InputIndexType, StreamlineType&);
  void RungeKuttaIntegrationNew(InputIndexType, StreamlineType&);
	
 
  void TrackingVolume();
	
// void PrintSelf(std::ostream& os, Indent indent) const;
 
private:
// This is only a temporary patch to achieve the program to compile in *** Windows:
#ifndef round
	int round( float val ){
		int r = (int)val;
		return val-r>0.5?r+1:r;
	}
#endif
	float m_CurvatureThreshold;
	float m_StepLength;
	float m_FaThreshold;
	StreamlineListType m_StreamlineList;
	unsigned int m_FibersPerVoxel;
	unsigned int m_MinimumFiberLength;
	CharImageTypePointer m_SeedsImage; 
	bool m_Connection;

};

}// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTractography.txx"

#endif

#endif

