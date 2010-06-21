// generated by Fast Light User Interface Designer (fluid) version 1.0107

#ifndef fltkImage2DViewerGUI_h
#define fltkImage2DViewerGUI_h
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <fltkImage2DViewerWindow.h>
#include <FL/Fl_Value_Slider.H>

class fltkImage2DViewerGUI {
public:
  fltkImage2DViewerGUI();
  Fl_Double_Window *externalWindow;
  fltk::Image2DViewerWindow *imageViewer;
  Fl_Double_Window *intensityWindow;
  Fl_Value_Slider *minimumSlider;
private:
  void cb_minimumSlider_i(Fl_Value_Slider*, void*);
  static void cb_minimumSlider(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *maximumSlider;
private:
  void cb_maximumSlider_i(Fl_Value_Slider*, void*);
  static void cb_maximumSlider(Fl_Value_Slider*, void*);
public:
  virtual ~fltkImage2DViewerGUI();
  void SetLabel(const char *label);
  void Show(void);
  void Hide(void);
  void Redraw(void);
  virtual void Update(void);
  virtual void RenderImage( double min, double max );
};
#endif
