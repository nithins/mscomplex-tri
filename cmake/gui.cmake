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

add_subdirectory(utls)

option(VIEWER_RENDER_AWESOME "build renderer to render awesome" OFF)

set(dmsc_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh.h

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/tri_edge.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/tri_edge.cpp

)

#set(spin_SRCS
#  ${CMAKE_CURRENT_SOURCE_DIR}/spin/spin_image.h
#  ${CMAKE_CURRENT_SOURCE_DIR}/spin/spin_image.cpp
#)

set(viewer_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer_mainwindow.h
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer_mainwindow.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/PythonQtScriptingConsole.h
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/PythonQtScriptingConsole.cpp

  )

QT4_WRAP_CPP(viewer_MOC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer_mainwindow.h
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/PythonQtScriptingConsole.h
  OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED
  )

set(viewer_UIS
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/trimesh_viewer_mainwindow.ui
  )

QT4_WRAP_UI(viewer_UI_SRCS ${viewer_UIS})


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
  ${PYTHONQT_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  )

configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/config.h)

add_executable(${PROJECT_NAME}_gui
  ${dmsc_SRCS}
  ${spin_SRCS}
  ${viewer_SRCS}
  ${viewer_MOC_SRCS}
  ${viewer_UI_SRCS}
  )

target_link_libraries(
  ${PROJECT_NAME}_gui
  utls
  ${Boost_LIBRARIES}
  ${QT_LIBRARIES}
  ${QT_ADDITIONAL_LIBRARIES}
  ${QGLVIEWER_LIBRARY}
  ${PYTHONQT_LIBRARY}
  ${PYTHON_LIBRARIES}
  GL
  )

