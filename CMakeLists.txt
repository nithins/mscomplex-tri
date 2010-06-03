#Define the project we are working on
project(mscomplex-2d)

# Look for various software we need
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Boost COMPONENTS program_options thread serialization regex REQUIRED)

set(QT_MIN_VERSION "4.5.0")
set(QT_USE_QTGUI TRUE)
set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTXML TRUE)

find_package(Qt4 REQUIRED)

include(${QT_USE_FILE})

find_package(QGLViewer REQUIRED)

SET(QXT_USE_QXTCORE TRUE)
SET(QXT_USE_QXTGUI TRUE)

find_package(Qxt REQUIRED)

find_package(PythonInterp REQUIRED)

find_package(OpenCL REQUIRED)

add_subdirectory(utls)

set(dmsc_SRCS

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.h

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/bitonic_sort.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/bitonic_sort.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/prefix_scan.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/prefix_scan.cpp
  )

QT4_ADD_RESOURCES(dmsc_RC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/ocl/mscomplex.qrc
  )

set(viewer_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.h
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.cpp
  )

QT4_WRAP_CPP(viewer_MOC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.h
  )

set(viewer_UIS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.ui
  )

QT4_WRAP_UI(viewer_UI_SRCS ${viewer_UIS})

file(GLOB viewer_shader_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/geom_1mfold.glsl
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/vert_1mfold.glsl
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/geom_2mfold.glsl
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/vert_2mfold.glsl
  )

set(viewer_shader_SRCS
  ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources.h
  ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources.cpp
  )

add_custom_command(
  OUTPUT ${viewer_shader_SRCS}
  COMMAND ${PYTHON_EXECUTABLE} ${utls_SOURCE_DIR}/copy_files_to_c_header.py
  ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources ${viewer_shader_FILES}
  DEPENDS ${viewer_shader_FILES} ${utls_SOURCE_DIR}/copy_files_to_c_header.py
  )

set(main_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/
  ${CMAKE_CURRENT_BINARY_DIR}/viewer/
  ${utls_SOURCE_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR}
  ${OPENCL_INCLUDE_DIRS}
  ${QT_ADDITIONAL_INCLUDE_PATHS}
  ${QGLVIEWER_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${QXT_INCLUDE_DIRS}
  )

add_executable(${PROJECT_NAME}
  ${dmsc_SRCS}
  ${dmsc_RC_SRCS}
  ${viewer_SRCS}
  ${viewer_MOC_SRCS}
  ${viewer_UI_SRCS}
  ${viewer_shader_SRCS}
  ${main_SRCS}
  )

target_link_libraries(
  ${PROJECT_NAME}
  utls
  ${Boost_LIBRARIES}
  ${OPENCL_LIBRARIES}
  ${QT_LIBRARIES}
  ${QT_ADDITIONAL_LIBRARIES}
  ${QGLVIEWER_LIBRARY}
  ${QXT_LIBRARIES}
  )
