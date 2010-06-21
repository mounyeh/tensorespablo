// generated by Fast Light User Interface Designer (fluid) version 1.0107

#ifndef fltkImageViewerGUI_h
#define fltkImageViewerGUI_h
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Choice.H>
#include <GLSliceView.h>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Browser.H>

class fltkImageViewerGUI {
public:
  fltkImageViewerGUI();
  virtual ~fltkImageViewerGUI();
  Fl_Double_Window* CreateGUI();
  Fl_Double_Window *iviewWindow;
  Fl_Group *glWindowGroup;
  Fl_Choice *orientationChoice;
private:
  void cb_orientationChoice_i(Fl_Choice*, void*);
  static void cb_orientationChoice(Fl_Choice*, void*);
  static Fl_Menu_Item menu_orientationChoice[];
  static Fl_Menu_Item menu_[];
  void cb_Value_i(Fl_Menu_*, void*);
  static void cb_Value(Fl_Menu_*, void*);
  void cb_Log_i(Fl_Menu_*, void*);
  static void cb_Log(Fl_Menu_*, void*);
  void cb_Opacity_i(Fl_Menu_*, void*);
  static void cb_Opacity(Fl_Menu_*, void*);
public:
  Fl_Value_Slider *sliceNumberSlider;
private:
  void cb_sliceNumberSlider_i(Fl_Value_Slider*, void*);
  static void cb_sliceNumberSlider(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *intensityWindowingMinSlider;
private:
  void cb_intensityWindowingMinSlider_i(Fl_Value_Slider*, void*);
  static void cb_intensityWindowingMinSlider(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *intensityWindowingMaxSlider;
private:
  void cb_intensityWindowingMaxSlider_i(Fl_Value_Slider*, void*);
  static void cb_intensityWindowingMaxSlider(Fl_Value_Slider*, void*);
  void cb_Zoom_i(Fl_Button*, void*);
  static void cb_Zoom(Fl_Button*, void*);
  void cb_Zoom1_i(Fl_Button*, void*);
  static void cb_Zoom1(Fl_Button*, void*);
  void cb_U_i(Fl_Button*, void*);
  static void cb_U(Fl_Button*, void*);
  void cb_Reset_i(Fl_Button*, void*);
  static void cb_Reset(Fl_Button*, void*);
  void cb_R_i(Fl_Button*, void*);
  static void cb_R(Fl_Button*, void*);
  void cb_L_i(Fl_Button*, void*);
  static void cb_L(Fl_Button*, void*);
  void cb_D_i(Fl_Button*, void*);
  static void cb_D(Fl_Button*, void*);
  void cb_Points_i(Fl_Button*, void*);
  static void cb_Points(Fl_Button*, void*);
public:
  Fl_Double_Window *clickedPointsWindow;
  Fl_Browser *clickedPointsBrowser;
private:
  void cb_Update_i(Fl_Button*, void*);
  static void cb_Update(Fl_Button*, void*);
  void cb_Clear_i(Fl_Button*, void*);
  static void cb_Clear(Fl_Button*, void*);
  void cb_Close_i(Fl_Button*, void*);
  static void cb_Close(Fl_Button*, void*);
public:
  Fl_Double_Window *overlayOpacityControlWindow;
private:
  void cb_Close1_i(Fl_Button*, void*);
  static void cb_Close1(Fl_Button*, void*);
public:
  Fl_Value_Slider *overlayOpacitySlider;
private:
  void cb_overlayOpacitySlider_i(Fl_Value_Slider*, void*);
  static void cb_overlayOpacitySlider(Fl_Value_Slider*, void*);
public:
  virtual void CreateGLSliceView( Fl_Group *,Fl_Gl_Window * w );
  virtual void AddMenuBarOptions(void);
  virtual void AddFilterMenuOptions();
  virtual void AddFileMenuOptions();
  virtual void SetImage( itk::ImageBase<3> * img );
  virtual void Show(void);
  virtual void Hide(void);
  virtual void Update(void);
  virtual void Synchronize(void);
  virtual void ImageMode(itk::ImageModeType mode);
  virtual void SelectSlice(unsigned int);
  virtual void SetIntensityWindowingMin(float);
  virtual void SetIntensityWindowingMax(float);
  virtual void CenterWindow(void);
  virtual void ZoomIn(void);
  virtual void ZoomOut(void);
  virtual void ShiftUp(void);
  virtual void ShiftDown(void);
  virtual void ShiftLeft(void);
  virtual void ShiftRight(void);
  virtual void SetOrientation(void);
  virtual void ShowClickedPointsWindow(void);
  virtual void ClearClickedPoints(void);
  virtual void ShowOverlayOpacityControl(void);
  virtual void SetOverlayOpacity(float);
  virtual void UpdateClickedPoints();
};
#endif
