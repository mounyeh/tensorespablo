// generated by Fast Light User Interface Designer (fluid) version 1.0107

#ifndef fltkDisplayGlWindowGUI_h
#define fltkDisplayGlWindowGUI_h
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <fltkGlWindowInteractive.h>
#include <fltkDrawer.h>

class fltkDisplayGlWindowGUI {
public:
  fltkDisplayGlWindowGUI();
private:
  Fl_Double_Window *parentWindow;
  fltk::GlWindowInteractive *drawWindow;
public:
  void Show(void);
  void SetLabel(const char *newlabel);
  void Redraw(void);
  void Size(unsigned int nx, unsigned int ny);
  void Hide(void);
  void Update(void);
  int GetWidth(void);
  int GetHeight(void);
  void MakeCurrent(void);
  int IsVisible(void);
  fltk::GlWindow::RedrawCommandType * GetRedrawCommand(void);
  itk::Object * GetNotifier(void);
  fltk::GlWindowInteractive * GetGlWindow(void);
  virtual void SaveImage(const char * filename);
};
#endif
