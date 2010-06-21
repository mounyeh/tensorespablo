// generated by Fast Light User Interface Designer (fluid) version 1.0107

#include "fltkTimeProbesCollectorGUI.h"

void TimeProbesCollectorGUI::cb_Close_i(Fl_Button*, void*) {
  this->Hide();
}
void TimeProbesCollectorGUI::cb_Close(Fl_Button* o, void* v) {
  ((TimeProbesCollectorGUI*)(o->parent()->user_data()))->cb_Close_i(o,v);
}

void TimeProbesCollectorGUI::cb_Clear_i(Fl_Button*, void*) {
  Clear();
}
void TimeProbesCollectorGUI::cb_Clear(Fl_Button* o, void* v) {
  ((TimeProbesCollectorGUI*)(o->parent()->user_data()))->cb_Clear_i(o,v);
}

void TimeProbesCollectorGUI::cb_Report_i(Fl_Button*, void*) {
  Report();
}
void TimeProbesCollectorGUI::cb_Report(Fl_Button* o, void* v) {
  ((TimeProbesCollectorGUI*)(o->parent()->user_data()))->cb_Report_i(o,v);
}

TimeProbesCollectorGUI::TimeProbesCollectorGUI() {
  Fl_Double_Window* w;
  { Fl_Double_Window* o = controlPanel = new Fl_Double_Window(378, 307, "Time Probes");
    w = o;
    o->box(FL_UP_BOX);
    o->user_data((void*)(this));
    { Fl_Button* o = new Fl_Button(295, 264, 75, 25, "Close");
      o->callback((Fl_Callback*)cb_Close);
    }
    { Fl_Button* o = new Fl_Button(15, 265, 73, 24, "Clear");
      o->callback((Fl_Callback*)cb_Clear);
    }
    { Fl_Button* o = new Fl_Button(100, 265, 73, 24, "Report");
      o->callback((Fl_Callback*)cb_Report);
    }
    { Fl_Check_Button* o = continuousCheckButton = new Fl_Check_Button(185, 264, 105, 25, "Continuous");
      o->box(FL_ENGRAVED_FRAME);
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->selection_color((Fl_Color)2);
    }
    { Fl_Scroll* o = new Fl_Scroll(10, 35, 360, 185);
      o->box(FL_DOWN_FRAME);
      { Fl_Pack* o = probesPack = new Fl_Pack(15, 40, 350, 180);
        o->end();
        Fl_Group::current()->resizable(o);
      }
      o->end();
      Fl_Group::current()->resizable(o);
    }
    { Fl_Box* o = new Fl_Box(20, 10, 95, 25, "Probe Tag");
      o->align(FL_ALIGN_CENTER|FL_ALIGN_INSIDE);
    }
    { Fl_Box* o = new Fl_Box(160, 10, 60, 25, "Starts");
      o->align(FL_ALIGN_CENTER|FL_ALIGN_INSIDE);
    }
    { Fl_Box* o = new Fl_Box(215, 10, 50, 25, "Stops");
      o->align(FL_ALIGN_CENTER|FL_ALIGN_INSIDE);
    }
    { Fl_Box* o = new Fl_Box(270, 10, 90, 25, "Time");
      o->align(FL_ALIGN_CENTER|FL_ALIGN_INSIDE);
    }
    timeResolution = new Fl_Value_Output(223, 230, 76, 25, "Time Resolution in this System :");
    new Fl_Box(300, 230, 65, 25, "seconds");
    o->end();
  }
}

TimeProbesCollectorGUI::~TimeProbesCollectorGUI() {
}

void TimeProbesCollectorGUI::Show(void) {
  controlPanel->show();
}

void TimeProbesCollectorGUI::Hide(void) {
  controlPanel->hide();
}
