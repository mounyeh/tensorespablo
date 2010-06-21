#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkEventObject.h>

class GenericImageToImageFilter {

public: 
  
  typedef float                              InputPixelType;
  typedef itk::Image< InputPixelType, 3 >    InputImageType;
  typedef InputImageType::Pointer  InputImagePointerType;
  typedef itk::ImageToImageFilter<InputImageType,InputImageType> itkFilter;
  //typedef itk::MemberCommand<UsimagToolConsole> CommandType;
  //typedef itk::SmartPointer<CommandType> CommandPointer;
  GenericImageToImageFilter();
  ~GenericImageToImageFilter();
  void SetFilter(itkFilter*);
  void SetInput(InputImageType*);
  void SetParameters();
  void SetNumberofThreads(int);
  void Update(void);
  //void AddObserver(const EventObject &event,CommandPointer callback);
  void* GetOutput(void);
  itkFilter* myfilter;

};
