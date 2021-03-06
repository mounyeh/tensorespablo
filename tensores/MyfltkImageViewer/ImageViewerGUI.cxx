// generated by Fast Light User Interface Designer (fluid) version 1.0107

#include "ImageViewerGUI.h"

void ImageViewerGUI::cb__i(Fl_Button*, void*) {
  this->ZoomIn();
}
void ImageViewerGUI::cb_(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb__i(o,v);
}

void ImageViewerGUI::cb_1_i(Fl_Button*, void*) {
  this->ZoomOut();
}
void ImageViewerGUI::cb_1(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_1_i(o,v);
}

void ImageViewerGUI::cb_2_i(Fl_Button*, void*) {
  ShiftUp();
}
void ImageViewerGUI::cb_2(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_2_i(o,v);
}

void ImageViewerGUI::cb_center_i(Fl_Button*, void*) {
  CenterWindow();
}
void ImageViewerGUI::cb_center(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_center_i(o,v);
}

void ImageViewerGUI::cb_3_i(Fl_Button*, void*) {
  ShiftRight();
}
void ImageViewerGUI::cb_3(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_3_i(o,v);
}

void ImageViewerGUI::cb_4_i(Fl_Button*, void*) {
  ShiftLeft();
}
void ImageViewerGUI::cb_4(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_4_i(o,v);
}

void ImageViewerGUI::cb_21_i(Fl_Button*, void*) {
  ShiftDown();
}
void ImageViewerGUI::cb_21(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_21_i(o,v);
}

void ImageViewerGUI::cb_sliceNumberSlider_i(Fl_Value_Slider* o, void*) {
  this->SelectSlice(o->value());
}
void ImageViewerGUI::cb_sliceNumberSlider(Fl_Value_Slider* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_sliceNumberSlider_i(o,v);
}

void ImageViewerGUI::cb_T_i(Fl_Button*, void*) {
  //this->transpose(); 
Update();
}
void ImageViewerGUI::cb_T(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_T_i(o,v);
}

void ImageViewerGUI::cb_v_i(Fl_Button*, void*) {
  bool t = this->viewValue();
this->ViewValue(!t); 
Update();
}
void ImageViewerGUI::cb_v(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_v_i(o,v);
}

void ImageViewerGUI::cb_d_i(Fl_Button*, void*) {
  bool t = this->viewDetails();
this->ViewDetails(!t); 
Update();
}
void ImageViewerGUI::cb_d(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_d_i(o,v);
}

void ImageViewerGUI::cb_c_i(Fl_Button*, void*) {
  bool t = this->viewCrosshairs();
this->ViewCrosshairs(!t); 
Update();
}
void ImageViewerGUI::cb_c(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_c_i(o,v);
}

void ImageViewerGUI::cb_intensityWindowingMaxSlider_i(Fl_Value_Slider* o, void*) {
  SetIntensityWindowingMax(o->value());
}
void ImageViewerGUI::cb_intensityWindowingMaxSlider(Fl_Value_Slider* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_intensityWindowingMaxSlider_i(o,v);
}

void ImageViewerGUI::cb_intensityWindowingMinSlider_i(Fl_Value_Slider* o, void*) {
  SetIntensityWindowingMin(o->value());
}
void ImageViewerGUI::cb_intensityWindowingMinSlider(Fl_Value_Slider* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_intensityWindowingMinSlider_i(o,v);
}

void ImageViewerGUI::cb_orientationChoice_i(Fl_Choice*, void*) {
  SetOrientation();
}
void ImageViewerGUI::cb_orientationChoice(Fl_Choice* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_orientationChoice_i(o,v);
}

Fl_Menu_Item ImageViewerGUI::menu_orientationChoice[] = {
 {"X", 0,  0, 0, 0, FL_NORMAL_LABEL, 0, 14, 0},
 {"Y", 0,  0, 0, 0, FL_NORMAL_LABEL, 0, 14, 0},
 {"Z", 0,  0, 0, 0, FL_NORMAL_LABEL, 0, 14, 0},
 {0,0,0,0,0,0,0,0,0}
};

void ImageViewerGUI::cb_x_i(Fl_Button*, void*) {
  bool t = this->flipX();
this->flipX(!t); 
Update();
}
void ImageViewerGUI::cb_x(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_x_i(o,v);
}

void ImageViewerGUI::cb_y_i(Fl_Button*, void*) {
  bool t = this->flipY();
this->flipY(!t); 
Update();
}
void ImageViewerGUI::cb_y(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_y_i(o,v);
}

void ImageViewerGUI::cb_z_i(Fl_Button*, void*) {
  bool t = this->flipZ();
this->flipZ(!t); 
Update();
}
void ImageViewerGUI::cb_z(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->parent()->user_data()))->cb_z_i(o,v);
}

void ImageViewerGUI::cb_Update_i(Fl_Button*, void*) {
  UpdateClickedPoints();
}
void ImageViewerGUI::cb_Update(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->user_data()))->cb_Update_i(o,v);
}

void ImageViewerGUI::cb_Clear_i(Fl_Button*, void*) {
  ClearClickedPoints();
}
void ImageViewerGUI::cb_Clear(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->user_data()))->cb_Clear_i(o,v);
}

void ImageViewerGUI::cb_Close_i(Fl_Button*, void*) {
  clickedPointsWindow->hide();
}
void ImageViewerGUI::cb_Close(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->user_data()))->cb_Close_i(o,v);
}

void ImageViewerGUI::cb_Close1_i(Fl_Button*, void*) {
  overlayOpacityControlWindow->hide();
}
void ImageViewerGUI::cb_Close1(Fl_Button* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->user_data()))->cb_Close1_i(o,v);
}

void ImageViewerGUI::cb_overlayOpacitySlider_i(Fl_Value_Slider* o, void*) {
  SetOverlayOpacity( o->value() );
}
void ImageViewerGUI::cb_overlayOpacitySlider(Fl_Value_Slider* o, void* v) {
  ((ImageViewerGUI*)(o->parent()->user_data()))->cb_overlayOpacitySlider_i(o,v);
}

ImageViewerGUI::ImageViewerGUI() {
}

ImageViewerGUI::~ImageViewerGUI() {
}

Fl_Double_Window* ImageViewerGUI::CreateGUI() {
  Fl_Double_Window* w;
  { Fl_Double_Window* o = iviewWindow = new Fl_Double_Window(260, 347);
    w = o;
    o->box(FL_PLASTIC_DOWN_BOX);
    o->color(FL_FOREGROUND_COLOR);
    o->labelsize(10);
    o->user_data((void*)(this));
    { Fl_Group* o = glWindowGroup = new Fl_Group(0, 0, 260, 345);
      { Fl_Gl_Window* o = new Fl_Gl_Window(0, 0, 260, 250, "3D Win");
        o->box(FL_PLASTIC_DOWN_BOX);
        o->color(FL_BACKGROUND_COLOR);
        o->selection_color(FL_BACKGROUND_COLOR);
        o->labeltype(FL_NORMAL_LABEL);
        o->labelfont(0);
        o->labelsize(14);
        o->labelcolor(FL_FOREGROUND_COLOR);
        o->align(FL_ALIGN_CENTER);
        o->when(FL_WHEN_RELEASE);
        CreateGLSliceView( glWindowGroup, o );
      }
      o->end();
    }
    { Fl_Group* o = new Fl_Group(0, 250, 260, 96);
      { Fl_Button* o = new Fl_Button(3, 254, 15, 15, "+");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_);
      }
      { Fl_Button* o = new Fl_Button(3, 269, 15, 15, "-");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_1);
      }
      { Fl_Button* o = new Fl_Button(36, 255, 15, 15, "@2<");
        o->box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_2);
      }
      { Fl_Button* o = new Fl_Button(68, 259, 35, 20, "center");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_center);
      }
      { Fl_Button* o = new Fl_Button(51, 263, 15, 15, "@>");
        o->box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_3);
      }
      { Fl_Button* o = new Fl_Button(21, 263, 15, 15, "@<");
        o->box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_4);
      }
      { Fl_Button* o = new Fl_Button(36, 270, 15, 15, "@2>");
        o->box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_21);
      }
      { Fl_Value_Slider* o = sliceNumberSlider = new Fl_Value_Slider(34, 289, 215, 15, "slice");
        o->type(1);
        o->box(FL_ENGRAVED_BOX);
        o->color(FL_INACTIVE_COLOR);
        o->selection_color(FL_SELECTION_COLOR);
        o->labelfont(1);
        o->labelsize(12);
        o->maximum(100);
        o->step(1);
        o->textsize(14);
        o->callback((Fl_Callback*)cb_sliceNumberSlider);
        o->align(132);
      }
      { Fl_Button* o = new Fl_Button(109, 254, 15, 15, "T");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelfont(1);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_T);
      }
      { Fl_Button* o = new Fl_Button(126, 254, 15, 15, "v");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelfont(1);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_v);
      }
      { Fl_Button* o = new Fl_Button(144, 254, 15, 15, "d");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelfont(1);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_d);
      }
      { Fl_Button* o = new Fl_Button(162, 254, 15, 15, "c");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelfont(1);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_c);
      }
      { Fl_Value_Slider* o = intensityWindowingMaxSlider = new Fl_Value_Slider(34, 308, 215, 15, "max");
        o->type(1);
        o->box(FL_ENGRAVED_BOX);
        o->color((Fl_Color)22);
        o->selection_color(FL_SELECTION_COLOR);
        o->labelfont(1);
        o->labelsize(12);
        o->textsize(14);
        o->callback((Fl_Callback*)cb_intensityWindowingMaxSlider);
        o->align(FL_ALIGN_LEFT);
      }
      { Fl_Value_Slider* o = intensityWindowingMinSlider = new Fl_Value_Slider(34, 325, 215, 15, "min");
        o->type(1);
        o->box(FL_ENGRAVED_BOX);
        o->color((Fl_Color)22);
        o->selection_color(FL_SELECTION_COLOR);
        o->labelfont(1);
        o->labelsize(12);
        o->textsize(14);
        o->callback((Fl_Callback*)cb_intensityWindowingMinSlider);
        o->align(FL_ALIGN_LEFT);
      }
      { Fl_Choice* o = orientationChoice = new Fl_Choice(198, 259, 42, 16);
        o->box(FL_PLASTIC_THIN_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color(FL_INACTIVE_COLOR);
        o->labelfont(1);
        o->labelsize(10);
        o->textsize(10);
        o->callback((Fl_Callback*)cb_orientationChoice);
        o->menu(menu_orientationChoice);
      }
      { Fl_Button* o = new Fl_Button(128, 272, 15, 15, "x");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelfont(1);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_x);
      }
      { Fl_Button* o = new Fl_Button(145, 272, 15, 15, "y");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelfont(1);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_y);
      }
      { Fl_Button* o = new Fl_Button(162, 272, 15, 15, "z");
        o->box(FL_PLASTIC_UP_BOX);
        o->down_box(FL_PLASTIC_DOWN_BOX);
        o->color((Fl_Color)4);
        o->selection_color((Fl_Color)4);
        o->labelfont(1);
        o->labelsize(10);
        o->callback((Fl_Callback*)cb_z);
      }
      { Fl_Text_Display* o = new Fl_Text_Display(127, 272, 51, 17, "Flip");
        o->box(FL_BORDER_FRAME);
        o->labelfont(1);
        o->labelsize(12);
        o->align(FL_ALIGN_LEFT);
      }
      o->end();
    }
    o->set_modal();
    o->end();
    o->resizable(o);
  }
  { Fl_Double_Window* o = clickedPointsWindow = new Fl_Double_Window(290, 405, "Clicked Points");
    w = o;
    o->user_data((void*)(this));
    { Fl_Browser* o = clickedPointsBrowser = new Fl_Browser(10, 15, 270, 345);
      o->textfont(4);
    }
    { Fl_Button* o = new Fl_Button(27, 370, 65, 25, "Update");
      o->callback((Fl_Callback*)cb_Update);
    }
    { Fl_Button* o = new Fl_Button(115, 370, 65, 25, "Clear");
      o->callback((Fl_Callback*)cb_Clear);
    }
    { Fl_Button* o = new Fl_Button(200, 370, 65, 25, "Close");
      o->callback((Fl_Callback*)cb_Close);
    }
    o->end();
  }
  { Fl_Double_Window* o = overlayOpacityControlWindow = new Fl_Double_Window(380, 70, "Overlay Opacity");
    w = o;
    o->box(FL_PLASTIC_UP_BOX);
    o->color((Fl_Color)59);
    o->labelsize(12);
    o->user_data((void*)(this));
    { Fl_Button* o = new Fl_Button(135, 44, 110, 16, "Close");
      o->box(FL_PLASTIC_UP_BOX);
      o->down_box(FL_PLASTIC_DOWN_BOX);
      o->color(FL_FOREGROUND_COLOR);
      o->labelsize(12);
      o->callback((Fl_Callback*)cb_Close1);
    }
    { Fl_Value_Slider* o = overlayOpacitySlider = new Fl_Value_Slider(10, 19, 360, 16, "Opacity");
      o->type(3);
      o->box(FL_PLASTIC_UP_BOX);
      o->color(FL_FOREGROUND_COLOR);
      o->selection_color((Fl_Color)21);
      o->labelsize(10);
      o->value(0.5);
      o->callback((Fl_Callback*)cb_overlayOpacitySlider);
      o->align(FL_ALIGN_TOP);
    }
    o->end();
  }
  return w;
}

void ImageViewerGUI::CreateGLSliceView( Fl_Group *,Fl_Gl_Window * w ) {
}

void ImageViewerGUI::AddMenuBarOptions(void) {
}

void ImageViewerGUI::AddFilterMenuOptions() {
}

void ImageViewerGUI::AddFileMenuOptions() {
}

void ImageViewerGUI::SetImage( itk::ImageBase<3> * img ) {
}

void ImageViewerGUI::Show(void) {
}

void ImageViewerGUI::Hide(void) {
}

void ImageViewerGUI::Update(void) {
}

void ImageViewerGUI::Synchronize(void) {
}

void ImageViewerGUI::ImageMode(itk::ImageModeType mode) {
}

void ImageViewerGUI::SelectSlice(unsigned int) {
}

void ImageViewerGUI::SetIntensityWindowingMin(float) {
}

void ImageViewerGUI::SetIntensityWindowingMax(float) {
}

void ImageViewerGUI::CenterWindow(void) {
}

void ImageViewerGUI::ZoomIn(void) {
}

void ImageViewerGUI::ZoomOut(void) {
}

void ImageViewerGUI::ShiftUp(void) {
}

void ImageViewerGUI::ShiftDown(void) {
}

void ImageViewerGUI::ShiftLeft(void) {
}

void ImageViewerGUI::ShiftRight(void) {
}

void ImageViewerGUI::SetOrientation(void) {
}

void ImageViewerGUI::ShowClickedPointsWindow(void) {
}

void ImageViewerGUI::ClearClickedPoints(void) {
}

void ImageViewerGUI::ShowOverlayOpacityControl(void) {
}

void ImageViewerGUI::SetOverlayOpacity(float) {
}

void ImageViewerGUI::UpdateClickedPoints() {
}

void ImageViewerGUI::ViewDetails(bool) {
}

void ImageViewerGUI::ViewValue(bool) {
}

void ImageViewerGUI::ViewCrosshairs(bool) {
}

bool ImageViewerGUI::viewCrosshairs(void) {
  return true;
}

bool ImageViewerGUI::viewValue(void) {
  return true;
}

bool ImageViewerGUI::viewDetails(void) {
  return true;
}

bool ImageViewerGUI::flipX(void) {
  return true;
}

void ImageViewerGUI::flipX(bool) {
}

bool ImageViewerGUI::flipY(void) {
  return true;
}

void ImageViewerGUI::flipY(bool) {
}

bool ImageViewerGUI::flipZ(void) {
  return true;
}

void ImageViewerGUI::flipZ(bool) {
}
