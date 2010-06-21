/*=========================================================================

  Program:   UsimagTool
  Language:  C++
  Date:      5-07-2007
  Version:   1.0

  Copyright (c) 2007 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/

#include "UsimagToolConsole.h"

int main()
{

 
  UsimagToolConsole *console = new UsimagToolConsole();

  console->Show();
  
  Fl::run();
  
  delete console;

  return 0;
}



