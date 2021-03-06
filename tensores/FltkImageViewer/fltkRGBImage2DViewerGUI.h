// generated by Fast Light User Interface Designer (fluid) version 1.0107

#ifndef fltkRGBImage2DViewerGUI_h
#define fltkRGBImage2DViewerGUI_h
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <fltkRGBImage2DViewerWindow.h>
#include <FL/Fl_Value_Slider.H>

class fltkRGBImage2DViewerGUI {
public:
  fltkRGBImage2DViewerGUI();
  Fl_Double_Window *externalWindow;
  fltk::RGBImage2DViewerWindow *imageViewer;
  Fl_Double_Window *intensityWindow;
private:
  void cb_Min_i(Fl_Value_Slider*, void*);
  static void cb_Min(Fl_Value_Slider*, void*);
  void cb_Max_i(Fl_Value_Slider*, void*);
  static void cb_Max(Fl_Value_Slider*, void*);
public:
  virtual ~fltkRGBImage2DViewerGUI();
  void SetLabel(const char *label);
  void Show(void);
  void Hide(void);
  void Redraw(void);
  virtual void SetMin(double val);
  virtual void SetMax(double val);
  virtual void Update(void);
};
#endif
