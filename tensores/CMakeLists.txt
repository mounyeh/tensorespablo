PROJECT( UsimagTool )

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build UsimagTool without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK REQUIRED)
IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
  SET(VTK_LIBRARIES
  vtkFLTK
  vtkRendering 
  vtkGraphics 
  vtkHybrid 
  vtkImaging 
  vtkIO 
  vtkFiltering 
  vtkCommon 
  vtkWidgets
  vtkVolumeRendering)
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build UsimagTool without VTK. Please set VTK_DIR.")
ENDIF(VTK_FOUND)

FIND_PACKAGE(FLTK REQUIRED)
IF(FLTK_FOUND)
  INCLUDE_DIRECTORIES(${FLTK_INCLUDE_DIR})
ELSE(FLTK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build UsimagTool without FLTK. Please set FLTK_DIR.")
ENDIF(FLTK_FOUND)

IF(WIN32)
   ADD_DEFINITIONS(-DWIN32)
   ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
   SET(GL_LIBRARIES)
ENDIF(WIN32)

IF(UNIX)
   ADD_DEFINITIONS(-DUNIX=1)
   SET(GL_LIBRARIES GL GLU)
ENDIF(UNIX)

# Macro for adding resources to GUI application on the Mac
#
#IF(APPLE AND NOT FLTK_USE_X)
#    FIND_PROGRAM(ITK_APPLE_RESOURCE Rez /Developer/Tools)
#   FIND_FILE(ITK_FLTK_RESOURCE mac.r 
#      ${FLTK_DIR}/FL)
#    IF(NOT ITK_FLTK_RESOURCE)
#      MESSAGE("Fltk resources not found, GUI application will not respond #to mouse events")
#    ENDIF(NOT ITK_FLTK_RESOURCE)
#
#    MACRO(ADD_GUI_EXECUTABLE name sources)
#      ADD_EXECUTABLE(${name} ${sources})
#      INSTALL_TARGETS(/bin ${name})
#      SET(EXEC_PATH ${EXECUTABLE_OUTPUT_PATH})
#      IF(NOT EXEC_PATH)
#        SET(EXEC_PATH ${CMAKE_CURRENT_BINARY_DIR})
#      ENDIF(NOT EXEC_PATH)
#        IF(ITK_APPLE_RESOURCE)
#          ADD_CUSTOM_COMMAND(SOURCE ${name}
#                             COMMAND ${ITK_APPLE_RESOURCE}
#                             ARGS -t APPL ${ITK_FLTK_RESOURCE} -o
#                             ${EXEC_PATH}/${name}
#                             TARGET ${name})
#        ENDIF(ITK_APPLE_RESOURCE)
#    ENDMACRO(ADD_GUI_EXECUTABLE)
#ELSE(APPLE AND NOT FLTK_USE_X)
#  MACRO(ADD_GUI_EXECUTABLE name sources)
#    ADD_EXECUTABLE(${name} ${sources})
#    INSTALL_TARGETS(/bin ${name})
#  ENDMACRO(ADD_GUI_EXECUTABLE)
#ENDIF(APPLE AND NOT FLTK_USE_X)

IF(APPLE)
 ADD_DEFINITIONS(-DAPPLE=1)

   FIND_PROGRAM(APPLE_RESOURCE Rez /Developer/Tools)
   IF(APPLE_RESOURCE)
     ADD_CUSTOM_COMMAND(TARGET UsimagTool POST_BUILD
                        COMMAND ${APPLE_RESOURCE} Carbon.r -o ./UsimagTool)
   ENDIF(APPLE_RESOURCE)
ENDIF(APPLE)

LINK_DIRECTORIES(
.
)

INCLUDE_DIRECTORIES(
.
MyfltkImageViewer
FltkImageViewer
VtkFltk
ASR
SRAD
Demons3D
vtkMarcacionElipse
wiener
knn-1canal
3DViewer
IO
)

#### FltkImageViewer 
SET(FltkImageViewerDIR FltkImageViewer)

SET(FltkImageViewer_SRCS
  ${FltkImageViewerDIR}/fltkBox3D.cxx
  ${FltkImageViewerDIR}/fltkCone3D.cxx
  ${FltkImageViewerDIR}/fltkClippingPlane3DDrawer.cxx
  ${FltkImageViewerDIR}/fltkCylinder3D.cxx
  ${FltkImageViewerDIR}/fltkFrustum3D.cxx
  ${FltkImageViewerDIR}/fltkGlWindow.cxx
  ${FltkImageViewerDIR}/fltkGlWindowInteractive.cxx
  ${FltkImageViewerDIR}/fltkGlWindowInteractiveFlat.cxx
  ${FltkImageViewerDIR}/fltkImage2DViewerWindow.cxx
  ${FltkImageViewerDIR}/fltkRGBImage2DViewerWindow.cxx
  ${FltkImageViewerDIR}/fltkLightButton.cxx
  ${FltkImageViewerDIR}/fltkProgressBar.cxx
  ${FltkImageViewerDIR}/fltkShape3D.cxx
  ${FltkImageViewerDIR}/fltkSphere3D.cxx
  ${FltkImageViewerDIR}/fltkWindowInteractive.cxx
  ${FltkImageViewerDIR}/fltkTimeProbeGUI.cxx
  ${FltkImageViewerDIR}/fltkTimeProbesCollector.cxx
  ${FltkImageViewerDIR}/fltkUtils.cxx
)

SET(FltkImageViewer_GUI_SRCS
  ${FltkImageViewerDIR}/fltkClippingPlane3DDrawerGUI.fl
  ${FltkImageViewerDIR}/fltkDisplayGlWindowGUI.fl
  ${FltkImageViewerDIR}/fltkDisplayGlWindowFlatGUI.fl
  ${FltkImageViewerDIR}/fltkFrustumFunctionControlGUI.fl
  ${FltkImageViewerDIR}/fltkImageViewerGUI.fl
  ${FltkImageViewerDIR}/fltkImage2DViewerGUI.fl
  ${FltkImageViewerDIR}/fltkRGBImage2DViewerGUI.fl
  ${FltkImageViewerDIR}/fltkRegularStepGradientDescentGUI.fl
  ${FltkImageViewerDIR}/fltkSlice2DDrawerGUI.fl
  ${FltkImageViewerDIR}/fltkSlice3DDrawerGUI.fl
  ${FltkImageViewerDIR}/fltkSphereFunctionControlGUI.fl
  ${FltkImageViewerDIR}/fltkTimeProbesCollectorGUI.fl
)


FLTK_WRAP_UI(ITKFltkImageViewer ${FltkImageViewer_GUI_SRCS} )
ADD_LIBRARY(ITKFltkImageViewer ${FltkImageViewer_SRCS} ${ITKFltkImageViewer_FLTK_UI_SRCS})
#### FltkImageViewer 

#### VTkFltkImageViewer 
FLTK_WRAP_UI(ITKVtkFltk VtkFltk/fltkVTKImageViewerGUI.fl)

ADD_LIBRARY(ITKVtkFltk
  VtkFltk/vtkFlRenderWindowInteractor.cxx
  VtkFltk/fltkVTKImageViewerBase.cxx
  ${ITKVtkFltk_FLTK_UI_SRCS}
)
#### VtkFltkImageViewer 

ADD_LIBRARY(vtkFLTK 
3DViewer/Fl_VTK_Window.cxx
3DViewer/vtkFLTKObjectFactory.cxx
3DViewer/vtkFLTKOpenGLRenderWindow.cxx
3DViewer/vtkFLTKRenderWindowInteractor.cxx
)

SET (Kit_LIBS
  ${VTK_FLTK_REQUIRED_KITS}
  ${FLTK_LIBRARIES}
  ${OPENGL_LIBRARIES}
)

FLTK_WRAP_UI(UsimagTool 
UsimagToolGUI.fl 
TensorGUI.fl
SegmentationGUI.fl
FilteringGUI.fl
BasicOpGUI.fl
MyfltkImageViewer/ImageViewerGUI.fl
)

SET(UsimagTool_SRCS
UsimagTool.cxx 
UsimagToolBase.cxx 
UsimagToolConsole.cxx 
TensorConsole.cxx
SegmentationConsole.cxx
FilteringConsole.cxx
BasicOpConsole.cxx
VolumesContainer.txx
geodesicPath3D.txx
GenericImageToImageFilter.cxx
MyfltkImageViewer/MyfltkImageViewer.txx 
MyfltkImageViewer/MyfltkColorImageViewer.txx 
MyfltkImageViewer/ImageViewer.txx 
Kretz/itkImageIOFactory.cxx 
Kretz/itkVOLImageIO.cxx 
Kretz/itkVOLImageIOFactory.cxx 
Kretz/itkTensorVTKImageIO.cxx 
Kretz/itkTensorVTKImageIOFactory.cxx 
Kretz/importvol.cxx 
vtkMarcacionElipse/vtkMarcacionElipsoide.cxx 
vtkMarcacionElipse/vtkPlantillaAjustada.cxx 
vtkMarcacionElipse/vtkFuncionVerosimilitud.cxx 
vtkMarcacionElipse/vtkOptimizaContorno.cxx 
vtkMarcacionElipse/vtkValidacion.cxx 
wiener/itkWienerFilter.txx
knn-1canal/classifyKNNCore.txx
knn-1canal/voronoiFilter.cxx
3DViewer/Viewer3D.cxx
vtkTensorGlyphDTI.cxx
strain/vtkTensorGlyphStrain.cxx
)

#Variable for choosing development or reduced interface
OPTION(PHYSICIANS_INTERFACE "Interface for clinicians" ON)
IF(PHYSICIANS_INTERFACE)
  SET_SOURCE_FILES_PROPERTIES(${UsimagTool_SRCS} COMPILE_FLAGS -DFORPHYSICIANS=1)
ENDIF(PHYSICIANS_INTERFACE)

ADD_EXECUTABLE(UsimagTool 
${UsimagTool_FLTK_UI_SRCS}
${UsimagTool_SRCS}
)



TARGET_LINK_LIBRARIES(UsimagTool 
${ITK_LIBRARIES}
ITKFltkImageViewer 
ITKVtkFltk 
${GL_LIBRARIES}
${FLTK_LIBRARIES}
${VTK_LIBRARIES}
)


