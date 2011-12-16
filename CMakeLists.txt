#Define the project we are working on
project(mscomplex-tri)

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

add_subdirectory(utls)

option(VIEWER_RENDER_AWESOME "build renderer to render awesome" OFF)

set(dmsc_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh.h

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_datamanager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_datamanager.cpp
  )

set(spin_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/spin/spin_image.h
  ${CMAKE_CURRENT_SOURCE_DIR}/spin/spin_image.cpp
)

set(viewer_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer_mainwindow.h
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer_mainwindow.cpp
  )

QT4_WRAP_CPP(viewer_MOC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer_mainwindow.h
  )

set(viewer_UIS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer_mainwindow.ui
  )

QT4_WRAP_UI(viewer_UI_SRCS ${viewer_UIS})

set(main_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/
  ${CMAKE_CURRENT_SOURCE_DIR}/spin/
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/
  ${CMAKE_CURRENT_BINARY_DIR}/viewer/
  ${utls_SOURCE_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/utls/
  ${QT_ADDITIONAL_INCLUDE_PATHS}
  ${QGLVIEWER_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${QXT_INCLUDE_DIRS}
  )

configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/config.h)

add_executable(${PROJECT_NAME}
  ${dmsc_SRCS}
  ${dmsc_RC_SRCS}
  ${spin_SRCS}
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
  ${QT_LIBRARIES}
  ${QT_ADDITIONAL_LIBRARIES}
  ${QGLVIEWER_LIBRARY}
  ${QXT_LIBRARIES}
  )
