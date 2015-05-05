cmake_minimum_required(VERSION 2.8)

project(HeightRidgeExtractionIn2D)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()

add_executable(HeightRidgeExtractor2DTest height_ridge_extractor_2d_test.cc
                                          height_ridge_extractor_2d.cc
                                          math.cc
                                          util.cc)
target_link_libraries(HeightRidgeExtractor2DTest ${VTK_LIBRARIES})

# add_executable(ScalarFieldTest scalar_field_test.cc
#                                util.cc)
# target_link_libraries(ScalarFieldTest ${VTK_LIBRARIES})

# add_executable(ITKTest itk_test.cc)
# target_link_libraries(ITKTest ${ITK_LIBRARIES})
# target_link_libraries(ITKTest ${VTK_LIBRARIES})

