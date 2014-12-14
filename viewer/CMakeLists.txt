project(mscomplex-tri-viewer)

cmake_minimum_required(VERSION 2.8.8)

set(QT_MIN_VERSION "4.5.0")
set(QT_USE_QTGUI TRUE)
set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTXML TRUE)

find_package(Qt4 REQUIRED)

include(${QT_USE_FILE})

find_package(QGLViewer REQUIRED)

find_package(PythonLibs REQUIRED)

find_package(PythonQt REQUIRED)

find_package(Eigen3 REQUIRED)

find_package(Boost COMPONENTS program_options serialization system REQUIRED)

option(VIEWER_RENDER_AWESOME "build renderer to render awesome" OFF)

#set(spin_SRCS
#  ${CMAKE_CURRENT_SOURCE_DIR}/spin/spin_image.h
#  ${CMAKE_CURRENT_SOURCE_DIR}/spin/spin_image.cpp
#)

add_subdirectory(utls)

set(viewer_SRCS
  trimesh_viewer.h
  main.cpp
  trimesh_viewer.cpp

  trimesh_viewer_mainwindow.h
  trimesh_viewer_mainwindow.cpp

  PythonQtScriptingConsole.h
  PythonQtScriptingConsole.cpp

  )

QT4_WRAP_CPP(viewer_MOC_SRCS
  trimesh_viewer_mainwindow.h
  PythonQtScriptingConsole.h
  OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED
  )

set(viewer_UIS
  trimesh_viewer_mainwindow.ui
  )

QT4_WRAP_UI(viewer_UI_SRCS ${viewer_UIS})


include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../core/
#  ${CMAKE_CURRENT_SOURCE_DIR}/../spin/ 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${utls_SOURCE_DIR}/include
  ${utls_BINARY_DIR}

  ${QT_ADDITIONAL_INCLUDE_PATHS}
  ${QGLVIEWER_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${PYTHONQT_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )

configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/config.h)

add_executable(mscomplex-tri-viewer
  $<TARGET_OBJECTS:mscomplex-tri-core>
#  ${spin_SRCS}
  ${viewer_SRCS}
  ${viewer_MOC_SRCS}
  ${viewer_UI_SRCS}
  )

target_link_libraries(
  mscomplex-tri-viewer
  utls
  ${Boost_LIBRARIES}
  ${QT_LIBRARIES}
  ${QT_ADDITIONAL_LIBRARIES}
  ${QGLVIEWER_LIBRARY}
  ${PYTHONQT_LIBRARY}
  ${PYTHON_LIBRARIES}
  GL
  )

