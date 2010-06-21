// generated by Fast Light User Interface Designer (fluid) version 1.0107

#include "fltkSlice3DDrawerGUI.h"

void fltkSlice3DDrawerGUI::cb_xScrollBar_i(Fl_Scrollbar*, void*) {
  SelectSliceX();
}
void fltkSlice3DDrawerGUI::cb_xScrollBar(Fl_Scrollbar* o, void* v) {
  ((fltkSlice3DDrawerGUI*)(o->parent()->user_data()))->cb_xScrollBar_i(o,v);
}

void fltkSlice3DDrawerGUI::cb_yScrollBar_i(Fl_Scrollbar*, void*) {
  SelectSliceY();
}
void fltkSlice3DDrawerGUI::cb_yScrollBar(Fl_Scrollbar* o, void* v) {
  ((fltkSlice3DDrawerGUI*)(o->parent()->user_data()))->cb_yScrollBar_i(o,v);
}

void fltkSlice3DDrawerGUI::cb_zScrollBar_i(Fl_Scrollbar*, void*) {
  SelectSliceZ();
}
void fltkSlice3DDrawerGUI::cb_zScrollBar(Fl_Scrollbar* o, void* v) {
  ((fltkSlice3DDrawerGUI*)(o->parent()->user_data()))->cb_zScrollBar_i(o,v);
}

void fltkSlice3DDrawerGUI::cb_xCheckButton_i(Fl_Check_Button*, void*) {
  SelectSliceX();
SelectSliceX();
}
void fltkSlice3DDrawerGUI::cb_xCheckButton(Fl_Check_Button* o, void* v) {
  ((fltkSlice3DDrawerGUI*)(o->parent()->user_data()))->cb_xCheckButton_i(o,v);
}

void fltkSlice3DDrawerGUI::cb_yCheckButton_i(Fl_Check_Button*, void*) {
  SelectSliceY();
SelectSliceY();
}
void fltkSlice3DDrawerGUI::cb_yCheckButton(Fl_Check_Button* o, void* v) {
  ((fltkSlice3DDrawerGUI*)(o->parent()->user_data()))->cb_yCheckButton_i(o,v);
}

void fltkSlice3DDrawerGUI::cb_zCheckButton_i(Fl_Check_Button*, void*) {
  SelectSliceZ();
SelectSliceZ();
}
void fltkSlice3DDrawerGUI::cb_zCheckButton(Fl_Check_Button* o, void* v) {
  ((fltkSlice3DDrawerGUI*)(o->parent()->user_data()))->cb_zCheckButton_i(o,v);
}

fltkSlice3DDrawerGUI::fltkSlice3DDrawerGUI() {
  Fl_Double_Window* w;
  { Fl_Double_Window* o = VolumeWindow = new Fl_Double_Window(452, 100, "Data Volume");
    w = o;
    o->user_data((void*)(this));
    { Fl_Scrollbar* o = xScrollBar = new Fl_Scrollbar(33, 22, 345, 17, "X : ");
      o->type(1);
      o->box(FL_DOWN_BOX);
      o->labelsize(10);
      o->callback((Fl_Callback*)cb_xScrollBar);
      o->align(FL_ALIGN_LEFT);
      o->when(FL_WHEN_RELEASE);
    }
    { Fl_Scrollbar* o = yScrollBar = new Fl_Scrollbar(33, 39, 345, 17, "Y : ");
      o->type(1);
      o->box(FL_DOWN_BOX);
      o->labelsize(12);
      o->callback((Fl_Callback*)cb_yScrollBar);
      o->align(FL_ALIGN_LEFT);
      o->when(FL_WHEN_RELEASE);
    }
    { Fl_Scrollbar* o = zScrollBar = new Fl_Scrollbar(33, 56, 345, 17, "Z : ");
      o->type(1);
      o->box(FL_DOWN_BOX);
      o->labelsize(10);
      o->callback((Fl_Callback*)cb_zScrollBar);
      o->align(FL_ALIGN_LEFT);
      o->when(FL_WHEN_RELEASE);
    }
    { Fl_Value_Output* o = xValueOutput = new Fl_Value_Output(379, 22, 34, 16);
      o->labelsize(10);
      o->textsize(10);
    }
    { Fl_Check_Button* o = xCheckButton = new Fl_Check_Button(412, 17, 30, 25);
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->selection_color((Fl_Color)2);
      o->callback((Fl_Callback*)cb_xCheckButton);
      o->align(FL_ALIGN_CENTER);
    }
    { Fl_Check_Button* o = yCheckButton = new Fl_Check_Button(412, 34, 26, 25);
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->selection_color((Fl_Color)2);
      o->callback((Fl_Callback*)cb_yCheckButton);
      o->align(FL_ALIGN_CENTER);
    }
    { Fl_Check_Button* o = zCheckButton = new Fl_Check_Button(412, 52, 24, 24);
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->selection_color((Fl_Color)2);
      o->callback((Fl_Callback*)cb_zCheckButton);
      o->align(FL_ALIGN_CENTER);
    }
    { Fl_Value_Output* o = yValueOutput = new Fl_Value_Output(379, 39, 34, 16);
      o->labelsize(10);
      o->textsize(10);
    }
    { Fl_Value_Output* o = zValueOutput = new Fl_Value_Output(379, 56, 34, 17);
      o->labelsize(10);
      o->textsize(10);
    }
    o->end();
  }
}

fltkSlice3DDrawerGUI::~fltkSlice3DDrawerGUI() {
}

void fltkSlice3DDrawerGUI::Show(void) {
}

void fltkSlice3DDrawerGUI::Hide(void) {
}

void fltkSlice3DDrawerGUI::SelectSliceX(void) {
}

void fltkSlice3DDrawerGUI::SelectSliceY(void) {
}

void fltkSlice3DDrawerGUI::SelectSliceZ(void) {
}

void fltkSlice3DDrawerGUI::SetLabel( const char * label ) {
  VolumeWindow->label( label );
}
