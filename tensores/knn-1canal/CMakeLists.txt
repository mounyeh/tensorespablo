# This project is designed to be built outside the Insight source tree.
PROJECT(ClassifyKNNGeodesic)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


ADD_EXECUTABLE(TestclassifyKNN TestclassifyKNN.cxx)
ADD_LIBRARY(voronoiFilter voronoiFilter.cxx)
ADD_LIBRARY(classifyKNNCore classifyKNNCore.cxx)
TARGET_LINK_LIBRARIES(TestclassifyKNN voronoiFilter classifyKNNCore ITKCommon ITKBasicFilters ITKIO )

