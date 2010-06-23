// generated by Fast Light User Interface Designer (fluid) version 1.0107

#ifndef fltkSphereFunctionControlGUI_h
#define fltkSphereFunctionControlGUI_h
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Adjuster.H>
#include <FL/Fl_Value_Output.H>
#include <FL/Fl_Box.H>

class fltkSphereFunctionControlGUI {
public:
  fltkSphereFunctionControlGUI();
  Fl_Double_Window *controlWindow;
  Fl_Adjuster *radiusAdjuster;
private:
  void cb_radiusAdjuster_i(Fl_Adjuster*, void*);
  static void cb_radiusAdjuster(Fl_Adjuster*, void*);
public:
  Fl_Value_Output *radiusValueOutput;
  Fl_Adjuster *xAdjuster;
private:
  void cb_xAdjuster_i(Fl_Adjuster*, void*);
  static void cb_xAdjuster(Fl_Adjuster*, void*);
public:
  Fl_Value_Output *xValueOutput;
  Fl_Adjuster *yAdjuster;
private:
  void cb_yAdjuster_i(Fl_Adjuster*, void*);
  static void cb_yAdjuster(Fl_Adjuster*, void*);
public:
  Fl_Value_Output *yValueOutput;
  Fl_Adjuster *zAdjuster;
private:
  void cb_zAdjuster_i(Fl_Adjuster*, void*);
  static void cb_zAdjuster(Fl_Adjuster*, void*);
public:
  Fl_Value_Output *zValueOutput;
  virtual ~fltkSphereFunctionControlGUI();
  virtual void SetRadius( double radius );
  virtual void SetCenterX( double x );
  virtual void SetCenterY( double y );
  virtual void SetCenterZ( double z );
  virtual void Show(void);
  virtual void Hide(void);
};
#endif