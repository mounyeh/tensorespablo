#include "GenericImageToImageFilter.h"

GenericImageToImageFilter::GenericImageToImageFilter() {
}

GenericImageToImageFilter::~GenericImageToImageFilter() {
}

void GenericImageToImageFilter::SetInput( InputImageType* input ) {
  this->myfilter->SetInput(input);
}

void GenericImageToImageFilter::SetFilter(itkFilter* filter) {
  this->myfilter = filter;
}

void GenericImageToImageFilter::SetNumberofThreads(int threads) {
  this->myfilter->SetNumberOfThreads(threads);
}

void GenericImageToImageFilter::SetParameters() {
}

void GenericImageToImageFilter::Update(void) {
  this->myfilter->Update();
}

void* GenericImageToImageFilter::GetOutput() {
  return this->myfilter->GetOutput();
}

//void* GenericImageToImageFilter::AddObserver(const EventObject &event,CommandPointer callback) {
//  return this->myfilter->AddObserver(event,callback);
//}
